#!/usr/bin/env python3
"""
Run Celligner on a union task matrix (genes × samples): split CPTAC/CCLE columns,
write temporary gene_matrix-style CSVs, call run_celligner_representation.

Usage:
  python scripts/benchmark/celligner_union_task.py \\
    --matrix data/processed_union/shared_gene_matrix_breast_subtype.csv \\
    --meta data/processed_union/sample_meta_breast_subtype.csv \\
    --out data/processed/methods/celligner/transformed_breast_subtype.csv \\
    --tmp-dir /tmp/celligner_union \\
    --no-umap
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "methods"))
from run_celligner_representation import run_celligner_representation  # noqa: E402


def load_union_genes_by_samples(path: Path) -> tuple[pd.DataFrame, str]:
    """Return genes × samples matrix; gene ID column name from first column."""
    df = pd.read_csv(path, header=0)
    gene_col = df.columns[0]
    df = df.set_index(gene_col)
    df.index = df.index.astype(str)
    return df, gene_col


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True)
    p.add_argument("--meta", required=True)
    p.add_argument("--out", required=True, help="Output genes × samples CSV path")
    p.add_argument("--tmp-dir", default=None, help="Temp workspace (default: system temp)")
    p.add_argument("--min-obs-frac", type=float, default=0.0)
    p.add_argument("--no-umap", action="store_true")
    args = p.parse_args()

    mat_path = Path(args.matrix)
    meta_path = Path(args.meta)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    mat, gene_col = load_union_genes_by_samples(mat_path)
    meta = pd.read_csv(meta_path)
    if not {"sample_id", "domain"}.issubset(meta.columns):
        raise SystemExit("meta needs sample_id, domain")

    cptac_ids = meta.loc[meta["domain"].str.upper() == "CPTAC", "sample_id"].astype(str)
    ccle_ids = meta.loc[meta["domain"].str.upper() == "CCLE", "sample_id"].astype(str)
    c_cols = [c for c in cptac_ids if c in mat.columns]
    e_cols = [c for c in ccle_ids if c in mat.columns]
    if len(c_cols) < 4 or len(e_cols) < 4:
        raise SystemExit(f"Too few domain columns CPTAC={len(c_cols)} CCLE={len(e_cols)}")

    tmp_root = Path(args.tmp_dir) if args.tmp_dir else Path(tempfile.mkdtemp(prefix="celligner_union_"))
    tmp_root.mkdir(parents=True, exist_ok=True)

    cptac_g = mat[c_cols].copy()
    ccle_g = mat[e_cols].copy()
    cptac_path = tmp_root / "cptac_genes_x_samples.csv"
    ccle_path = tmp_root / "ccle_genes_x_samples.csv"
    cptac_g.insert(0, gene_col, cptac_g.index)
    ccle_g.insert(0, gene_col, ccle_g.index)
    cptac_g.to_csv(cptac_path, index=False)
    ccle_g.to_csv(ccle_path, index=False)

    outdir = tmp_root / "celligner_out"
    res = run_celligner_representation(
        cptac_matrix_path=str(cptac_path),
        ccle_matrix_path=str(ccle_path),
        cptac_meta_path=str(meta_path),
        ccle_meta_path=None,
        outdir=str(outdir),
        min_obs_frac=args.min_obs_frac,
        impute_strategy="median",
        compute_umap=not args.no_umap,
    )

    src = Path(res["matrix_path"])
    shutil.copy2(src, out_path)
    print(f"[celligner_union_task] Wrote {out_path} (status={res.get('status')})")

    if args.tmp_dir is None:
        shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    main()
