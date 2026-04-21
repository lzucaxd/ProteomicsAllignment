#!/usr/bin/env python3
"""
Trial run: Celligner on a subsample of CPTAC breast + CCLE.

Takes a random but reproducible subsample of both domains,
runs Celligner, and evaluates basic metrics:
  - Did fit/transform succeed?
  - Domain effect before vs after (PC1 R²)
  - Marker direction preservation
  - Subtype separation within each domain

Usage:
  R_HOME=/Library/Frameworks/R.framework/Resources \
  R_PROFILE_USER=/dev/null \
  python scripts/benchmark/trial_celligner_subsample.py
"""
import os, sys, time
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "models" / "celligner-master"))

CPTAC_BREAST = REPO / "data" / "results" / "PDC000120" / "gene_matrix.csv"
CCLE_MATRIX  = REPO / "data" / "results" / "CCLE_corrected" / "gene_matrix.csv"
OUTDIR       = REPO / "reports" / "benchmark_master" / "celligner_trial"

SEED = 42
N_CPTAC = 100
N_CCLE  = 150
MIN_OBS_FRAC = 0.7

SUBTYPE_BASAL   = ["HCC70", "HCC1806", "HCC1143", "MDA-MB-468"]
SUBTYPE_LUMINAL = ["CAMA-1", "MCF7", "T-47D", "ZR-75-1"]

MARKERS = {
    "FOXA1": +1, "GATA3": +1, "ESR1": +1, "KRT5": -1,
    "KRT14": -1, "KRT17": -1, "EGFR": -1,
}


def load_matrix(path):
    df = pd.read_csv(path, index_col=0)
    if "UniProtID" in df.columns:
        df = df.drop(columns=["UniProtID"])
    return df.T.astype(float)


def compute_pc1_domain_r2(mat, domain_labels):
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(mat)
    domain_num = np.array([0 if d == "CPTAC" else 1 for d in domain_labels]).reshape(-1, 1)
    r2 = LinearRegression().fit(domain_num, pcs[:, 0]).score(domain_num, pcs[:, 0])
    return r2, pca.explained_variance_ratio_[:2]


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.RandomState(SEED)

    print("=" * 60)
    print("  CELLIGNER SUBSAMPLE TRIAL")
    print("=" * 60)

    # Load data
    print("\n[1/5] Loading matrices...")
    cptac_full = load_matrix(CPTAC_BREAST)
    ccle_full  = load_matrix(CCLE_MATRIX)
    print(f"  CPTAC: {cptac_full.shape[0]} samples × {cptac_full.shape[1]} genes")
    print(f"  CCLE:  {ccle_full.shape[0]} samples × {ccle_full.shape[1]} genes")

    # Subsample CPTAC
    cptac_idx = rng.choice(cptac_full.shape[0], size=min(N_CPTAC, cptac_full.shape[0]), replace=False)
    cptac_sub = cptac_full.iloc[cptac_idx]

    # For CCLE, always keep subtype lines + random others
    subtype_lines = SUBTYPE_BASAL + SUBTYPE_LUMINAL
    ccle_cols_lower = {c.replace("-", "").replace("_", "").upper(): c for c in ccle_full.index}
    kept_ccle = []
    for sl in subtype_lines:
        sl_norm = sl.replace("-", "").replace("_", "").upper()
        matches = [c for k, c in ccle_cols_lower.items() if sl_norm in k]
        kept_ccle.extend(matches[:1])
    remaining = [c for c in ccle_full.index if c not in kept_ccle]
    n_extra = min(N_CCLE - len(kept_ccle), len(remaining))
    extra = list(rng.choice(remaining, size=n_extra, replace=False))
    kept_ccle.extend(extra)
    ccle_sub = ccle_full.loc[kept_ccle]

    print(f"\n  Subsample — CPTAC: {cptac_sub.shape[0]}, CCLE: {ccle_sub.shape[0]}")

    # Intersect genes, filter, impute
    shared = sorted(set(cptac_sub.columns) & set(ccle_sub.columns))
    cptac_s = cptac_sub[shared]
    ccle_s  = ccle_sub[shared]
    obs_frac = (cptac_s.notna().mean(axis=0) >= MIN_OBS_FRAC) & \
               (ccle_s.notna().mean(axis=0) >= MIN_OBS_FRAC)
    genes_kept = [g for g in shared if obs_frac[g]]
    cptac_s = cptac_s[genes_kept].fillna(cptac_s[genes_kept].median(axis=0))
    ccle_s  = ccle_s[genes_kept].fillna(ccle_s[genes_kept].median(axis=0))
    print(f"  Shared genes after filter: {len(genes_kept)}")

    # Pre-alignment metrics
    print("\n[2/5] Pre-alignment metrics...")
    combined_pre = pd.concat([cptac_s, ccle_s])
    domains_pre = ["CPTAC"] * len(cptac_s) + ["CCLE"] * len(ccle_s)
    r2_pre, pve_pre = compute_pc1_domain_r2(combined_pre.values, domains_pre)
    print(f"  PC1 domain R² (pre):  {r2_pre:.4f}")
    print(f"  PC1 variance explained: {pve_pre[0]*100:.1f}%")

    # Run Celligner
    print("\n[3/5] Running Celligner...")
    t0 = time.time()
    from celligner import Celligner
    model = Celligner()

    print(f"  Fitting on CCLE ({ccle_s.shape[0]} × {ccle_s.shape[1]})...")
    model.fit(ccle_s)
    print(f"  Transforming CPTAC ({cptac_s.shape[0]} × {cptac_s.shape[1]})...")
    model.transform(cptac_s)
    elapsed = time.time() - t0
    print(f"  Celligner completed in {elapsed:.1f}s")

    combined_post = model.combined_output
    if combined_post is None:
        print("ERROR: combined_output is None")
        return

    print(f"  Combined output shape: {combined_post.shape}")

    # Post-alignment metrics
    print("\n[4/5] Post-alignment metrics...")
    cptac_ids = set(cptac_s.index)
    domains_post = ["CPTAC" if s in cptac_ids else "CCLE" for s in combined_post.index]
    r2_post, pve_post = compute_pc1_domain_r2(combined_post.values, domains_post)
    print(f"  PC1 domain R² (post): {r2_post:.4f}")
    print(f"  PC1 variance explained: {pve_post[0]*100:.1f}%")
    print(f"  Domain R² change: {r2_pre:.4f} → {r2_post:.4f} (Δ = {r2_post - r2_pre:+.4f})")

    # Marker check on CCLE subtype subset
    print("\n[5/5] Marker direction check (CCLE Basal vs Luminal)...")
    ccle_aligned = combined_post.loc[[s for s in combined_post.index if s not in cptac_ids]]
    basal_in_aligned = [s for s in ccle_aligned.index
                        if any(b.replace("-","").upper() in s.replace("-","").replace("_","").upper()
                               for b in SUBTYPE_BASAL)]
    luminal_in_aligned = [s for s in ccle_aligned.index
                          if any(l.replace("-","").upper() in s.replace("-","").replace("_","").upper()
                                 for l in SUBTYPE_LUMINAL)]

    marker_results = []
    for gene, expected_sign in MARKERS.items():
        if gene not in combined_post.columns:
            marker_results.append((gene, expected_sign, np.nan, "missing"))
            continue
        if len(basal_in_aligned) == 0 or len(luminal_in_aligned) == 0:
            marker_results.append((gene, expected_sign, np.nan, "no_subtype_samples"))
            continue
        lum_mean = ccle_aligned.loc[luminal_in_aligned, gene].mean()
        bas_mean = ccle_aligned.loc[basal_in_aligned, gene].mean()
        fc = lum_mean - bas_mean
        observed_sign = 1 if fc > 0 else -1
        correct = "OK" if observed_sign == expected_sign else "WRONG"
        marker_results.append((gene, expected_sign, fc, correct))
        print(f"  {gene}: FC={fc:+.3f} expected={'+' if expected_sign>0 else '-'} → {correct}")

    # Save results
    results_lines = [
        "CELLIGNER SUBSAMPLE TRIAL RESULTS",
        "=" * 40,
        f"Date: {pd.Timestamp.now()}",
        f"CPTAC subsample: {N_CPTAC} (from PDC000120)",
        f"CCLE subsample: {N_CCLE}",
        f"Genes (after filter): {len(genes_kept)}",
        f"Imputation: median per gene",
        f"Min observation fraction: {MIN_OBS_FRAC}",
        f"Random seed: {SEED}",
        "",
        "DOMAIN EFFECT",
        f"  PC1 domain R² (pre-alignment):  {r2_pre:.4f}",
        f"  PC1 domain R² (post-alignment): {r2_post:.4f}",
        f"  Change: {r2_post - r2_pre:+.4f}",
        f"  Interpretation: {'Domain effect reduced' if r2_post < r2_pre else 'Domain effect NOT reduced'}",
        "",
        f"Runtime: {elapsed:.1f}s",
        "",
        "MARKER DIRECTION CHECK (CCLE Luminal-Basal)",
    ]
    n_correct = sum(1 for _, _, _, s in marker_results if s == "OK")
    n_checked = sum(1 for _, _, _, s in marker_results if s in ("OK", "WRONG"))
    for gene, exp, fc, status in marker_results:
        results_lines.append(f"  {gene}: FC={fc:+.3f} expected={'+' if exp>0 else '-'} → {status}" if not np.isnan(fc) else f"  {gene}: {status}")
    results_lines.append(f"\n  Direction accuracy: {n_correct}/{n_checked}")

    results_path = OUTDIR / "trial_results.txt"
    with open(results_path, "w") as f:
        f.write("\n".join(results_lines))

    # Save aligned matrix
    combined_post.to_csv(OUTDIR / "trial_aligned_matrix.csv")

    # Save PCA data for plotting
    pca = PCA(n_components=2)
    pcs_post = pca.fit_transform(combined_post.values)
    pc_df = pd.DataFrame({
        "PC1": pcs_post[:, 0], "PC2": pcs_post[:, 1],
        "domain": domains_post, "sample": combined_post.index
    })
    pc_df.to_csv(OUTDIR / "trial_pca_post_alignment.csv", index=False)

    pcs_pre = PCA(n_components=2).fit_transform(combined_pre.values)
    pc_pre_df = pd.DataFrame({
        "PC1": pcs_pre[:, 0], "PC2": pcs_pre[:, 1],
        "domain": domains_pre, "sample": combined_pre.index
    })
    pc_pre_df.to_csv(OUTDIR / "trial_pca_pre_alignment.csv", index=False)

    marker_df = pd.DataFrame(marker_results, columns=["gene", "expected_sign", "FC", "status"])
    marker_df.to_csv(OUTDIR / "trial_marker_check.csv", index=False)

    print(f"\n{'=' * 60}")
    print(f"  TRIAL COMPLETE")
    print(f"  Outputs: {OUTDIR}")
    print(f"  Domain R²: {r2_pre:.4f} → {r2_post:.4f}")
    print(f"  Markers correct: {n_correct}/{n_checked}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
