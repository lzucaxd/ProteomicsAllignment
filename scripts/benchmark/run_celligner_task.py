#!/usr/bin/env python3
"""
Task-specific Celligner alignment on the task union matrix.

Goal: keep a constant (task-defined) gene set across methods for fair comparison.

Inputs (from preprocessing):
  data/processed/union/shared_gene_matrix_{task}.csv   (genes × samples; first col GeneSymbol)
  data/processed/union/sample_meta_{task}.csv         (sample_id, domain, condition, ...)

Outputs:
  data/processed/methods/celligner/transformed_{task}.csv   (genes × samples; GeneSymbol first col)
  reports/benchmark_master/celligner_task/{task}/celligner_aligned_matrix.csv  (samples × genes)
  reports/benchmark_master/celligner_task/{task}/run_summary.json

Notes:
  - Unlike run_celligner_all_data.py, this does NOT do an additional prevalence/SD filter;
    it trusts the task union construction to define the gene set.
  - Imputation: within-domain gene-median (only for remaining NAs).
  - Standardization: z-score per gene within each domain, then Celligner fit on CCLE, transform CPTAC.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


def read_union_matrix(path: Path) -> pd.DataFrame:
    dt = pd.read_csv(path)
    gene_col = dt.columns[0]
    mat = dt.set_index(gene_col)
    # Expect genes × samples; coerce numeric
    mat = mat.apply(pd.to_numeric, errors="coerce")
    mat.index = mat.index.astype(str)
    return mat


def zscore_by_gene(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=0)
    sd = df.std(axis=0)
    sd[sd == 0] = 1.0
    return (df - mu) / sd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", default=".", help="Repository root")
    ap.add_argument("--task", required=True, choices=["breast_subtype", "breast_vs_lung"])
    args = ap.parse_args()

    repo = Path(args.repo_root).resolve()
    task = args.task
    t0 = time.time()

    # Ensure Celligner import path
    import sys

    sys.path.insert(0, str(repo / "models" / "celligner-master"))
    from celligner import Celligner  # type: ignore

    union_dir = repo / "data" / "processed" / "union"
    mat_path = union_dir / f"shared_gene_matrix_{task}.csv"
    meta_path = union_dir / f"sample_meta_{task}.csv"

    if not mat_path.exists():
        raise FileNotFoundError(mat_path)
    if not meta_path.exists():
        raise FileNotFoundError(meta_path)

    mat_gxs = read_union_matrix(mat_path)  # genes × samples
    meta = pd.read_csv(meta_path)
    # normalize domain casing
    meta["domain"] = meta["domain"].astype(str)
    meta["sample_id"] = meta["sample_id"].astype(str)

    # Subset to columns present in matrix
    meta = meta[meta["sample_id"].isin(mat_gxs.columns)].copy()
    if meta.empty:
        raise ValueError("No sample_meta rows matched union matrix columns")

    # Split by domain (samples × genes for Celligner)
    cptac_ids = meta.loc[meta["domain"].str.upper() == "CPTAC", "sample_id"].tolist()
    ccle_ids = meta.loc[meta["domain"].str.upper() == "CCLE", "sample_id"].tolist()
    if len(cptac_ids) < 4 or len(ccle_ids) < 4:
        raise ValueError(f"Too few samples: CPTAC={len(cptac_ids)}, CCLE={len(ccle_ids)}")

    cptac = mat_gxs[cptac_ids].T  # samples × genes
    ccle = mat_gxs[ccle_ids].T

    # Impute within-domain gene median (over samples)
    na_cptac = int(cptac.isna().sum().sum())
    na_ccle = int(ccle.isna().sum().sum())
    cptac = cptac.fillna(cptac.median(axis=0))
    ccle = ccle.fillna(ccle.median(axis=0))

    # Z-score per gene, per domain
    cptac_z = zscore_by_gene(cptac)
    ccle_z = zscore_by_gene(ccle)

    # Run Celligner
    # The upstream Celligner defaults (pca_ncomp=70) assume large n; for task-specific
    # panels (e.g. CCLE subtype ~24 lines) we must reduce PCA dimensions.
    pca_ncomp = max(2, min(20, ccle_z.shape[0] - 1, ccle_z.shape[1] - 1))
    cpca_ncomp = max(2, min(4, pca_ncomp - 1))
    model = Celligner(pca_ncomp=pca_ncomp, cpca_ncomp=cpca_ncomp)
    model.fit(ccle_z)
    model.transform(cptac_z)
    combined_post = model.combined_output  # samples × genes

    # Write outputs
    out_task = repo / "reports" / "benchmark_master" / "celligner_task" / task
    out_task.mkdir(parents=True, exist_ok=True)
    combined_post.to_csv(out_task / "celligner_aligned_matrix.csv")

    # Also write transformed matrix into methods/ for the benchmark wrapper
    out_methods = repo / "data" / "processed" / "methods" / "celligner"
    out_methods.mkdir(parents=True, exist_ok=True)
    gxs_out = combined_post.T.copy()
    gxs_out.index.name = "GeneSymbol"
    gxs_out.to_csv(out_methods / f"transformed_{task}.csv")

    summary = {
        "task": task,
        "n_genes": int(mat_gxs.shape[0]),
        "n_samples_total": int(meta.shape[0]),
        "n_cptac": int(len(cptac_ids)),
        "n_ccle": int(len(ccle_ids)),
        "na_imputed_cptac": na_cptac,
        "na_imputed_ccle": na_ccle,
        "runtime_seconds": round(time.time() - t0, 2),
    }
    with (out_task / "run_summary.json").open("w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

