#!/usr/bin/env python3
"""
Method 2 — Celligner-Aligned Representation Wrapper

Wraps the local Celligner Python package (models/celligner-master/celligner)
to produce a benchmark-ready aligned matrix from CPTAC and CCLE gene matrices.

Celligner expects:
  - pandas DataFrames: samples (rows) × genes (columns)
  - Gene IDs: Ensembl (ENSG...) or gene symbols (the Python port works with
    whatever column names are provided as long as they overlap)
  - Values: log2(X+1) expression (designed for TPM; we adapt TMT abundances)
  - No NaN values in the input

Our adaptation:
  1. Load CPTAC and CCLE gene_matrix.csv (genes × samples → transpose to samples × genes)
  2. Inner-join on shared genes, drop genes with any NA, or impute per-gene median
  3. Fit Celligner on CCLE (reference = cell lines), transform CPTAC (target = tumors)
  4. Extract combined_output, save as benchmark-format CSV
  5. Save QC outputs (DE genes, cPCA loadings, UMAP if computed)

Usage:
  python3 scripts/methods/run_celligner_representation.py \
    --cptac_matrix data/results/PDC000120/gene_matrix.csv \
    --ccle_matrix  data/results/CCLE_corrected/gene_matrix.csv \
    --cptac_meta   <cptac_sample_meta.csv> \
    --ccle_meta    <ccle_sample_meta.csv> \
    --outdir       reports/benchmark_master/methods/celligner_style

Dependencies:
  - PYTHONPATH includes models/celligner-master (or pip install the mnnpy submodule).
  - R on PATH with Bioconductor limma (Celligner DE step calls Rscript; see celligner/limma.py).
  - pip install ".[celligner]" for scanpy, umap-learn, python-igraph, etc.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
CELLIGNER_ROOT = REPO_ROOT / "models" / "celligner-master"

sys.path.insert(0, str(CELLIGNER_ROOT))

def load_gene_matrix(path: str | Path) -> pd.DataFrame:
    """Load gene_matrix.csv (genes × samples) → samples × genes DataFrame."""
    df = pd.read_csv(path, index_col=0)
    if "UniProtID" in df.columns:
        df = df.drop(columns=["UniProtID"])
    return df.T.astype(float)


def prepare_inputs(
    cptac_df: pd.DataFrame,
    ccle_df: pd.DataFrame,
    min_obs_frac: float = 0.5,
    impute_strategy: str = "median",
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    """
    Intersect gene features, handle NAs, return Celligner-ready DataFrames.
    Returns (ccle_clean, cptac_clean, included_genes, excluded_genes).
    """
    shared_genes = sorted(set(cptac_df.columns) & set(ccle_df.columns))
    cptac_sub = cptac_df[shared_genes]
    ccle_sub = ccle_df[shared_genes]

    obs_cptac = cptac_sub.notna().mean(axis=0)
    obs_ccle = ccle_sub.notna().mean(axis=0)
    keep_mask = (obs_cptac >= min_obs_frac) & (obs_ccle >= min_obs_frac)

    included = [g for g in shared_genes if keep_mask[g]]
    excluded = [g for g in shared_genes if not keep_mask[g]]

    cptac_filt = cptac_sub[included]
    ccle_filt = ccle_sub[included]

    if impute_strategy == "median":
        cptac_filt = cptac_filt.fillna(cptac_filt.median(axis=0))
        ccle_filt = ccle_filt.fillna(ccle_filt.median(axis=0))
    elif impute_strategy == "drop":
        na_genes = cptac_filt.columns[cptac_filt.isna().any() | ccle_filt.isna().any()]
        excluded.extend(na_genes.tolist())
        included = [g for g in included if g not in set(na_genes)]
        cptac_filt = cptac_filt[included]
        ccle_filt = ccle_filt[included]

    return ccle_filt, cptac_filt, included, excluded


def run_celligner_representation(
    cptac_matrix_path: str,
    ccle_matrix_path: str,
    cptac_meta_path: str | None,
    ccle_meta_path: str | None,
    outdir: str,
    min_obs_frac: float = 0.5,
    impute_strategy: str = "median",
    compute_umap: bool = True,
) -> dict:
    """
    Run Celligner alignment and save benchmark-format outputs.

    Returns a dict mirroring the R method interface contract.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load matrices
    cptac_df = load_gene_matrix(cptac_matrix_path)
    ccle_df = load_gene_matrix(ccle_matrix_path)

    # Prepare inputs
    ccle_clean, cptac_clean, included_genes, excluded_genes = prepare_inputs(
        cptac_df, ccle_df, min_obs_frac=min_obs_frac, impute_strategy=impute_strategy
    )

    n_cptac = len(cptac_clean)
    n_ccle = len(ccle_clean)
    n_genes = len(included_genes)

    notes = [
        "Method: celligner",
        f"Date: {datetime.now().isoformat()}",
        "",
        "Description:",
        "  Celligner (Broad Institute) aligns tumor and cell line expression profiles",
        "  using contrastive PCA to remove domain-specific variation followed by mutual",
        "  nearest neighbors (MNN) batch correction. Originally designed for RNA-seq",
        "  log2(TPM+1); adapted here for TMT proteomics log2 abundance.",
        "",
        f"CPTAC samples: {n_cptac}",
        f"CCLE samples: {n_ccle}",
        f"Genes included: {n_genes}",
        f"Genes excluded: {len(excluded_genes)}",
        f"Imputation strategy: {impute_strategy}",
        "",
        f"Local Celligner repo: {CELLIGNER_ROOT}",
    ]

    qc_paths = {}

    try:
        from celligner import Celligner

        # Default PCA_NCOMP=70 exceeds n_samples for small CCLE arms (e.g. 25 lines);
        # scanpy PCA requires n_comps < min(n_samples, n_genes).
        n_min_domain = min(n_cptac, n_ccle)
        pca_ncomp = int(min(70, max(2, n_min_domain - 1), max(2, n_genes - 1)))
        use_low_mem = (n_genes > 6000) or (n_min_domain < 40)
        from celligner.params import TOP_K_GENES

        top_k_genes = min(TOP_K_GENES, n_genes)

        model = Celligner(
            topKGenes=top_k_genes,
            pca_ncomp=pca_ncomp,
            low_mem=use_low_mem,
        )

        notes.append("")
        notes.append("Celligner parameters (benchmark-adapted):")
        notes.append(f"  topKGenes: {model.topKGenes}")
        notes.append(f"  pca_ncomp: {model.pca_ncomp} (capped from default 70 for n_min_domain={n_min_domain})")
        notes.append(f"  low_mem: {model.low_mem}")
        notes.append(f"  cpca_ncomp: {model.cpca_ncomp}")
        notes.append(f"  mnn_kwargs: {model.mnn_kwargs}")

        print(f"[celligner] Fitting on CCLE ({n_ccle} samples × {n_genes} genes)...")
        model.fit(ccle_clean)

        print(f"[celligner] Transforming CPTAC ({n_cptac} samples × {n_genes} genes)...")
        model.transform(cptac_clean)

        combined = model.combined_output
        if combined is None:
            raise RuntimeError("Celligner transform() produced no combined_output")

        # Save DE genes
        if model.de_genes is not None:
            de_path = outdir / "celligner_de_genes.txt"
            with open(de_path, "w") as f:
                f.write("\n".join(model.de_genes))
            qc_paths["de_genes"] = str(de_path)

        # Save cPCA loadings
        if model.cpca_loadings is not None:
            cpca_df = pd.DataFrame(
                model.cpca_loadings,
                columns=included_genes[:model.cpca_loadings.shape[1]]
                if model.cpca_loadings.shape[1] <= n_genes
                else [f"gene_{i}" for i in range(model.cpca_loadings.shape[1])],
            )
            cpca_path = outdir / "celligner_cpca_loadings.csv"
            cpca_df.to_csv(cpca_path, index=False)
            qc_paths["cpca_loadings"] = str(cpca_path)

        # UMAP
        if compute_umap:
            try:
                print("[celligner] Computing UMAP...")
                model.computeMetricsForOutput(UMAP_only=True)
                if model.umap_reduced is not None:
                    umap_path = outdir / "celligner_umap.csv"
                    model.umap_reduced.to_csv(umap_path)
                    qc_paths["umap"] = str(umap_path)
            except Exception as e:
                warnings.warn(f"UMAP computation failed: {e}")
                notes.append(f"WARNING: UMAP computation failed: {e}")

        # Transpose back to genes × samples for benchmark contract
        matrix_out = combined.T
        matrix_out.index.name = "Gene"

        notes.extend([
            "",
            "Output:",
            f"  Combined matrix shape: {matrix_out.shape[0]} genes × {matrix_out.shape[1]} samples",
            "  The output is in Celligner-aligned space (cPCA-regressed + MNN-corrected).",
            "  This is NOT a TMT-native abundance matrix.",
            "  Absolute values are not directly comparable to MSstatsTMT log2 abundances.",
            "  Relative differences (fold-changes) within the aligned space are meaningful",
            "  for cross-domain comparison.",
            "",
            "What is preserved:",
            "  - Relative gene-level differences between samples within and across domains",
            "  - Broad cancer-type structure and tissue-of-origin signal",
            "",
            "What may be distorted:",
            "  - Absolute protein abundance levels",
            "  - Domain-specific biology that is confounded with batch (e.g., if a real",
            "    tissue vs cell line difference looks like a batch effect to cPCA)",
            "  - Subtype contrasts if subtype signal correlates with domain variation",
        ])

        status = "FULLY_IMPLEMENTED"

    except ImportError as e:
        notes.extend([
            "",
            f"ERROR: Could not import Celligner (or a dependency raised ImportError mid-run): {e}",
            "Python: pip install -e models/celligner-master/mnnpy && pip install 'harmonize[celligner]'",
            "R: BiocManager::install('limma') — DE uses Rscript (no rpy2).",
            "",
            "SCAFFOLD: Output matrix is a PLACEHOLDER (raw concatenation, not aligned).",
        ])
        status = "SCAFFOLDED_IMPORT_ERROR"

        combined_raw = pd.concat([cptac_clean, ccle_clean])
        matrix_out = combined_raw.T
        matrix_out.index.name = "Gene"

    except Exception as e:
        print(f"[celligner] Runtime error (scaffolding): {e!r}", file=sys.stderr)
        notes.extend([
            "",
            f"ERROR: Celligner execution failed: {e}",
            "SCAFFOLD: Output matrix is a PLACEHOLDER (raw concatenation, not aligned).",
        ])
        status = "SCAFFOLDED_RUNTIME_ERROR"

        combined_raw = pd.concat([cptac_clean, ccle_clean])
        matrix_out = combined_raw.T
        matrix_out.index.name = "Gene"

    # Save transformed matrix
    matrix_path = outdir / "transformed_matrix.csv"
    matrix_out.to_csv(matrix_path)

    # Build sample metadata
    sample_ids = list(matrix_out.columns)
    cptac_ids = set(cptac_clean.index)
    sample_meta = pd.DataFrame({
        "sample_id": sample_ids,
        "domain": ["CPTAC" if s in cptac_ids else "CCLE" for s in sample_ids],
        "condition": "Sample",
    })
    if cptac_meta_path and os.path.exists(cptac_meta_path):
        ext_meta = pd.read_csv(cptac_meta_path)
        sample_meta = sample_meta.merge(ext_meta, on="sample_id", how="left", suffixes=("", "_ext"))
    if ccle_meta_path and os.path.exists(ccle_meta_path):
        ext_meta = pd.read_csv(ccle_meta_path)
        sample_meta = sample_meta.merge(ext_meta, on="sample_id", how="left", suffixes=("", "_ext"))
    sample_meta.to_csv(outdir / "sample_metadata.csv", index=False)

    # Feature metadata
    feature_meta = pd.DataFrame({
        "gene": included_genes + excluded_genes,
        "included": [True] * len(included_genes) + [False] * len(excluded_genes),
        "exclusion_reason": [None] * len(included_genes) + ["low_obs_or_na"] * len(excluded_genes),
    })
    feature_meta.to_csv(outdir / "feature_metadata.csv", index=False)

    # Save notes
    notes_path = outdir / "method_notes.txt"
    with open(notes_path, "w") as f:
        f.write("\n".join(notes))
    qc_paths["notes"] = str(notes_path)

    # Save QC paths
    with open(outdir / "qc_paths.json", "w") as f:
        json.dump(qc_paths, f, indent=2)

    print(f"[celligner] Status: {status}")
    print(f"[celligner] Outputs saved to {outdir}")

    return {
        "matrix_path": str(matrix_path),
        "sample_meta_path": str(outdir / "sample_metadata.csv"),
        "feature_meta_path": str(outdir / "feature_metadata.csv"),
        "method_name": "celligner",
        "status": status,
        "notes": notes,
        "qc_paths": qc_paths,
    }


def main():
    parser = argparse.ArgumentParser(description="Method 2: Celligner-aligned representation")
    parser.add_argument("--cptac_matrix", required=True, help="Path to CPTAC gene_matrix.csv")
    parser.add_argument("--ccle_matrix", required=True, help="Path to CCLE gene_matrix.csv")
    parser.add_argument("--cptac_meta", default=None, help="Path to CPTAC sample metadata CSV")
    parser.add_argument("--ccle_meta", default=None, help="Path to CCLE sample metadata CSV")
    parser.add_argument("--outdir", default="reports/benchmark_master/methods/celligner_style")
    parser.add_argument("--min_obs_frac", type=float, default=0.5)
    parser.add_argument("--impute", default="median", choices=["median", "drop"])
    parser.add_argument("--no_umap", action="store_true")
    args = parser.parse_args()

    run_celligner_representation(
        cptac_matrix_path=args.cptac_matrix,
        ccle_matrix_path=args.ccle_matrix,
        cptac_meta_path=args.cptac_meta,
        ccle_meta_path=args.ccle_meta,
        outdir=args.outdir,
        min_obs_frac=args.min_obs_frac,
        impute_strategy=args.impute,
        compute_umap=not args.no_umap,
    )


if __name__ == "__main__":
    main()
