#!/usr/bin/env python3
"""
Sanity-check script for CCLE MSstatsTMT output (gene matrix, annotation, QC).
Validates matrix shape, missing values, zero-variance genes, sample intensity,
PCA, and annotation consistency before downstream analyses.

Usage:
  python data/scripts/check_ccle_matrix.py
  python data/scripts/check_ccle_matrix.py --matrix data/results/CCLE/gene_matrix.csv

Use the same Python that has pandas/numpy/matplotlib/sklearn installed (e.g. python3.12 if
packages are installed there and `python3` points elsewhere).

Outputs:
  data/results/CCLE/pca_plot.png  (PC1 vs PC2 scatter)
"""
import argparse
import os
import sys

try:
    import numpy as np
    import pandas as pd
except ImportError as e:
    print("Missing dependency:", e, file=sys.stderr)
    print("Install with: pip install pandas numpy matplotlib scikit-learn", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Sanity-check CCLE gene matrix and annotation")
    parser.add_argument(
        "--matrix",
        default="",
        help="Path to gene_matrix.csv (default: data/results/CCLE/gene_matrix.csv)",
    )
    parser.add_argument("--annotation", default="", help="Path to annotation_filled.csv")
    parser.add_argument("--outdir", default="", help="Directory for pca_plot.png")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_ccle = os.path.normpath(os.path.join(script_dir, "..", "results", "CCLE"))
    matrix_path = os.path.abspath(args.matrix) if args.matrix and os.path.isfile(args.matrix) else (os.path.abspath(args.matrix) if args.matrix else os.path.join(data_ccle, "gene_matrix.csv"))
    if not os.path.isfile(matrix_path) and args.matrix:
        matrix_path = os.path.join(data_ccle, "gene_matrix.csv")
    ann_path = os.path.abspath(args.annotation) if args.annotation and os.path.isfile(args.annotation) else os.path.join(data_ccle, "annotation_filled.csv")
    outdir = os.path.abspath(args.outdir) if args.outdir and os.path.isdir(args.outdir) else data_ccle
    if args.annotation and not os.path.isfile(ann_path):
        ann_path = os.path.abspath(args.annotation)
    if args.outdir and not os.path.isdir(outdir):
        outdir = os.path.abspath(args.outdir)

    if not os.path.isfile(matrix_path):
        print(f"ERROR: Matrix not found: {matrix_path}")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # 1. Load the gene matrix
    # -------------------------------------------------------------------------
    df = pd.read_csv(matrix_path)
    gene_col = df.columns[0]
    # Gene matrix may have GeneSymbol, UniProtID, then sample columns
    if len(df.columns) > 1 and df.columns[1] == "UniProtID":
        sample_cols = df.columns[2:].tolist()
    else:
        sample_cols = df.columns[1:].tolist()
    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    n_genes, n_samples = mat.shape

    print("--- 1. Matrix shape ---")
    print("Matrix shape:", mat.shape)
    print("Number of genes:", n_genes)
    print("Number of samples:", n_samples)

    if n_genes < 8000 or n_genes > 12000:
        print(f"  WARNING: Genes ({n_genes}) outside expected range 8000–12000")
    if n_samples < 300 or n_samples > 400:
        print(f"  WARNING: Samples ({n_samples}) outside expected range 300–400")

    # -------------------------------------------------------------------------
    # 2. Missing value statistics
    # -------------------------------------------------------------------------
    overall_missing = mat.isna().sum().sum() / (n_genes * n_samples) * 100
    missing_per_gene = mat.isna().mean(axis=1) * 100
    missing_per_sample = mat.isna().mean(axis=0) * 100

    print("\n--- 2. Missing value statistics ---")
    print(f"Overall missing %: {overall_missing:.2f}%")
    print(f"Median missing per gene: {missing_per_gene.median():.2f}%")
    print(f"Median missing per sample: {missing_per_sample.median():.2f}%")

    if overall_missing > 60:
        print(f"  WARNING: Overall missing ({overall_missing:.1f}%) > 60%")

    # -------------------------------------------------------------------------
    # 3. Zero / constant genes
    # -------------------------------------------------------------------------
    with np.errstate(invalid="ignore"):
        gene_var = mat.var(axis=1, skipna=True)
    zero_var_genes = (gene_var == 0).sum()
    all_na_genes = mat.isna().all(axis=1).sum()

    print("\n--- 3. Zero / constant genes ---")
    print("Genes with zero variance:", int(zero_var_genes))
    print("Genes with all NA:", int(all_na_genes))

    # -------------------------------------------------------------------------
    # 4. Sample intensity distribution
    # -------------------------------------------------------------------------
    # Values are already log2 abundance from MSstatsTMT
    median_per_sample = mat.median(axis=0, skipna=True)
    var_per_sample = mat.var(axis=0, skipna=True)
    sample_median_mean = median_per_sample.mean()
    sample_median_std = median_per_sample.std()
    if pd.isna(sample_median_std) or sample_median_std == 0:
        sample_median_std = 1e-10
    outlier_mask = (median_per_sample < sample_median_mean - 3 * sample_median_std) | (
        median_per_sample > sample_median_mean + 3 * sample_median_std
    )
    n_outlier_samples = outlier_mask.sum()

    print("\n--- 4. Sample intensity distribution (log2) ---")
    print("Median intensity per sample: mean = {:.3f}, std = {:.3f}".format(sample_median_mean, sample_median_std))
    print("Variance per sample: mean = {:.3f}".format(var_per_sample.mean()))
    if n_outlier_samples > 0:
        print(f"  FLAG: {int(n_outlier_samples)} sample(s) with median outside 3 SD: {median_per_sample[outlier_mask].index.tolist()[:5]}{'...' if n_outlier_samples > 5 else ''}")

    # -------------------------------------------------------------------------
    # 5. PCA sanity check
    # -------------------------------------------------------------------------
    var_explained = [0.0] * 5
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        mat_centered = mat - mat.mean(axis=0)
        mat_filled = mat_centered.fillna(0)
        scaler = StandardScaler(with_mean=False)
        X = scaler.fit_transform(mat_filled.T)
        pca = PCA(n_components=5)
        pca.fit(X)
        var_explained = (pca.explained_variance_ratio_ * 100).tolist()

        print("\n--- 5. PCA (n_components=5) ---")
        for i, v in enumerate(var_explained, 1):
            print(f"  PC{i} variance explained: {v:.2f}%")

        # Save PC1 vs PC2 scatter
        os.makedirs(outdir, exist_ok=True)
        pca_plot_path = os.path.join(outdir, "pca_plot.png")
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            scores = pca.transform(X)
            plt.figure(figsize=(7, 6))
            plt.scatter(scores[:, 0], scores[:, 1], alpha=0.6, s=20)
            plt.xlabel(f"PC1 ({var_explained[0]:.1f}%)")
            plt.ylabel(f"PC2 ({var_explained[1]:.1f}%)")
            plt.title("CCLE gene matrix: PC1 vs PC2")
            plt.tight_layout()
            plt.savefig(pca_plot_path, dpi=150)
            plt.close()
            print(f"  PCA plot saved: {pca_plot_path}")
        except Exception as e:
            print(f"  Could not save PCA plot: {e}")
    except ImportError as e:
        print("\n--- 5. PCA ---")
        print("  sklearn not available:", e)

    # -------------------------------------------------------------------------
    # 6. Annotation consistency
    # -------------------------------------------------------------------------
    if os.path.isfile(ann_path):
        ann = pd.read_csv(ann_path)
        if "BioReplicate" in ann.columns:
            # Only samples (exclude Norm/bridge) for comparison
            ann_samples = set(ann.loc[ann["Condition"].str.lower() != "norm", "BioReplicate"].dropna().unique())
        else:
            ann_samples = set()
        matrix_samples = set(sample_cols)
        in_ann_not_matrix = ann_samples - matrix_samples
        in_matrix_not_ann = matrix_samples - ann_samples

        print("\n--- 6. Annotation consistency ---")
        print("Unique samples in annotation (non-Norm):", len(ann_samples))
        print("Unique samples in gene matrix:", len(matrix_samples))
        if in_ann_not_matrix or in_matrix_not_ann:
            print("  WARNING: Sample name mismatch")
            if in_ann_not_matrix:
                print(f"    In annotation but not in matrix: {len(in_ann_not_matrix)} (e.g. {list(in_ann_not_matrix)[:3]})")
            if in_matrix_not_ann:
                print(f"    In matrix but not in annotation: {len(in_matrix_not_ann)} (e.g. {list(in_matrix_not_ann)[:3]})")
        else:
            print("  Sample names match.")
    else:
        print("\n--- 6. Annotation consistency ---")
        print(f"  Annotation file not found: {ann_path}")

    # -------------------------------------------------------------------------
    # 7. Bridge channels removed
    # -------------------------------------------------------------------------
    bridge_keywords = ["norm", "bridge", "pool"]
    found_bridge = [c for c in sample_cols if any(k in (c or "").lower() for k in bridge_keywords)]
    if found_bridge:
        print("\n--- 7. Bridge channels ---")
        print(f"  WARNING: Bridge-like names detected in sample names: {found_bridge}. Bridge channels may not have been removed.")
    else:
        print("\n--- 7. Bridge channels ---")
        print("  Bridge channels detected: NO (no Norm/Bridge/POOL in sample names)")

    # -------------------------------------------------------------------------
    # 8. Final summary report
    # -------------------------------------------------------------------------
    print("\n" + "=" * 50)
    print("CCLE PROTEOMICS SANITY CHECK")
    print("=" * 50)
    print(f"Genes: {n_genes}")
    print(f"Samples: {n_samples}")
    print(f"Missing %: {overall_missing:.2f}%")
    print(f"Zero variance genes: {int(zero_var_genes)}")
    print(f"Samples with outlier intensity: {int(n_outlier_samples)}")
    print("Bridge channels detected: NO" if not found_bridge else "Bridge channels detected: YES (WARN)")
    if len(var_explained) >= 5 and any(v > 0 for v in var_explained):
        print("PCA variance explained:")
        for i, v in enumerate(var_explained[:5], 1):
            print(f"  PC{i}: {v:.2f}%")
    print("=" * 50)


if __name__ == "__main__":
    main()
