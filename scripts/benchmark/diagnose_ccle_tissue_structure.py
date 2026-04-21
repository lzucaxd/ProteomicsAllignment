#!/usr/bin/env python3
"""
Diagnostic: Does CCLE proteomics have tissue-of-origin signal at all?
- PCA and UMAP on CCLE alone, colored by tissue
- Variance explained by tissue (ANOVA R² on PC1-5)
- Compare to CPTAC tissue structure
"""
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
OUTDIR = REPO / "reports/benchmark_master/diagnostics/ccle_tissue_check"
OUTDIR.mkdir(parents=True, exist_ok=True)

SITE_COLORS = {
    "Breast": "#E91E63", "Lung": "#4CAF50", "Ovary": "#00BCD4",
    "Uterine": "#F06292", "Colorectal": "#FF9800",
    "Blood/Lymphoid": "#9C27B0", "Skin": "#795548", "CNS": "#3F51B5",
    "Pancreas": "#CDDC39", "Upper GI": "#FF5722", "Endometrium": "#AB47BC",
    "Liver": "#8D6E63", "Kidney": "#5C6BC0", "Urinary Tract": "#26A69A",
    "Head & Neck": "#78909C", "Prostate": "#AED581", "Other": "#BDBDBD",
    "Unknown": "#E0E0E0",
}
TISSUE_REMAP = {
    "Breast": "Breast", "Lung": "Lung", "Ovary": "Ovary",
    "Large Intestine": "Colorectal",
    "Haematopoietic and Lymphoid Tissue": "Blood/Lymphoid",
    "Skin": "Skin", "Central Nervous System": "CNS",
    "Pancreas": "Pancreas", "Stomach": "Upper GI", "Oesophagus": "Upper GI",
    "Endometrium": "Endometrium", "Liver": "Liver", "Kidney": "Kidney",
    "Urinary Tract": "Urinary Tract", "Upper Aerodigestive Tract": "Head & Neck",
    "Prostate": "Prostate", "Soft Tissue": "Other", "Bone": "Other",
    "Thyroid": "Other", "Pleura": "Other", "Autonomic Ganglia": "Other",
    "Biliary Tract": "Other",
}

def load_gm(path):
    df = pd.read_csv(path, index_col=0)
    if "UniProtID" in df.columns:
        df = df.drop(columns=["UniProtID"])
    return df.T.astype(float)

def get_ccle_tissue_map():
    info = pd.read_csv(REPO / "data/ccle_peptide/sample_info_ccle.csv")
    raw = dict(zip(info["Cell Line"], info["Tissue of Origin"]))
    norm = {k.replace("-","").replace("_","").upper(): v for k, v in raw.items()}
    return raw, norm

def remap(t):
    return TISSUE_REMAP.get(t, "Other")

def tissue_r2(pcs, tissues):
    """R² from one-hot tissue encoding predicting each PC."""
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    t_enc = le.fit_transform(tissues)
    dummies = np.eye(len(le.classes_))[t_enc]
    r2s = []
    for j in range(pcs.shape[1]):
        lr = LinearRegression().fit(dummies, pcs[:, j])
        r2s.append(lr.score(dummies, pcs[:, j]))
    return r2s

def make_scatter(coords, tissues, pve, title, outpath, is_umap=False):
    fig, ax = plt.subplots(figsize=(12, 9))
    from collections import Counter
    ct = Counter(tissues)
    present = sorted(ct.keys(), key=lambda s: -ct[s])
    for site in present:
        idx = [i for i, t in enumerate(tissues) if t == site]
        c = SITE_COLORS.get(site, "gray")
        ax.scatter(coords[idx, 0], coords[idx, 1], c=c, alpha=0.7, s=40,
                   edgecolors="k", linewidths=0.3, label=f"{site} (n={ct[site]})")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, title="Tissue")
    if not is_umap:
        ax.set_xlabel(f"PC1 ({pve[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pve[1]*100:.1f}%)")
    else:
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved: {outpath.name}")


def main():
    print("=" * 70)
    print("  CCLE TISSUE STRUCTURE DIAGNOSTIC")
    print("=" * 70)

    # ── Load CCLE ────────────────────────────────────────────────────────
    ccle = load_gm(REPO / "data/results/CCLE_corrected/gene_matrix.csv")
    raw_map, norm_map = get_ccle_tissue_map()

    tissues = []
    for s in ccle.index:
        if s in raw_map:
            tissues.append(remap(raw_map[s]))
        else:
            s_n = s.replace("-","").replace("_","").upper()
            tissues.append(remap(norm_map.get(s_n, "Unknown")))
    tissues = np.array(tissues)

    # Prevalence + impute + z-score (same as alignment pipeline)
    prev = ccle.notna().mean(axis=0)
    keep = prev >= 0.70
    ccle_f = ccle.loc[:, keep]
    sd = ccle_f.std(axis=0, skipna=True)
    ccle_f = ccle_f.loc[:, sd >= 0.01]
    print(f"  CCLE: {ccle_f.shape[0]} samples × {ccle_f.shape[1]} genes (after 70% prev + SD filter)")

    na_before = ccle_f.isna().sum().sum()
    ccle_f = ccle_f.fillna(ccle_f.median(axis=0))
    print(f"  NAs imputed: {na_before}")

    # ── Analysis 1: RAW (no z-score) ────────────────────────────────────
    print("\n── CCLE raw (no z-score) ──")
    pca_raw = PCA(n_components=10)
    pcs_raw = pca_raw.fit_transform(ccle_f.values)
    pve_raw = pca_raw.explained_variance_ratio_
    r2_raw = tissue_r2(pcs_raw, tissues)
    print(f"  PVE: {[f'{v*100:.1f}%' for v in pve_raw[:5]]}")
    print(f"  Tissue R² per PC: {[f'{v:.4f}' for v in r2_raw[:5]]}")
    print(f"  Total tissue R² (PC1-5): {sum(r2_raw[:5] * pve_raw[:5]) / sum(pve_raw[:5]):.4f}")

    make_scatter(pcs_raw, tissues, pve_raw, "CCLE Raw PCA (no z-score)", OUTDIR / "ccle_pca_raw.png")

    # ── Analysis 2: Z-scored ────────────────────────────────────────────
    print("\n── CCLE z-scored ──")
    mu = ccle_f.mean(axis=0)
    s = ccle_f.std(axis=0); s[s == 0] = 1
    ccle_z = (ccle_f - mu) / s

    pca_z = PCA(n_components=10)
    pcs_z = pca_z.fit_transform(ccle_z.values)
    pve_z = pca_z.explained_variance_ratio_
    r2_z = tissue_r2(pcs_z, tissues)
    print(f"  PVE: {[f'{v*100:.1f}%' for v in pve_z[:5]]}")
    print(f"  Tissue R² per PC: {[f'{v:.4f}' for v in r2_z[:5]]}")
    print(f"  Total tissue R² (PC1-5): {sum(r2_z[:5] * pve_z[:5]) / sum(pve_z[:5]):.4f}")

    make_scatter(pcs_z, tissues, pve_z, "CCLE Z-scored PCA", OUTDIR / "ccle_pca_zscore.png")

    # ── Analysis 3: UMAP on CCLE alone ──────────────────────────────────
    print("\n── CCLE UMAP ──")
    import umap
    um = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=42, n_jobs=1)
    umap_coords = um.fit_transform(ccle_z.values)
    make_scatter(umap_coords, tissues, None, "CCLE UMAP (z-scored)", OUTDIR / "ccle_umap.png", is_umap=True)

    # ── Compare: CPTAC tissue structure ──────────────────────────────────
    print("\n── CPTAC tissue structure (for comparison) ──")
    studies = {
        "PDC000120": "Breast", "PDC000127": "Ovary",
        "PDC000153": "Lung", "PDC000204": "Uterine",
    }
    parts = []
    cptac_tissues = []
    for sid, cancer in studies.items():
        gm = load_gm(REPO / f"data/results/{sid}/gene_matrix.csv")
        parts.append(gm)
        cptac_tissues.extend([cancer] * len(gm))
    cptac = pd.concat(parts, join="outer")
    cptac_tissues = np.array(cptac_tissues)

    prev_c = cptac.notna().mean(axis=0)
    cptac_f = cptac.loc[:, prev_c >= 0.70]
    sd_c = cptac_f.std(axis=0, skipna=True)
    cptac_f = cptac_f.loc[:, sd_c >= 0.01]
    cptac_f = cptac_f.fillna(cptac_f.median(axis=0))
    mu_c = cptac_f.mean(axis=0)
    s_c = cptac_f.std(axis=0); s_c[s_c == 0] = 1
    cptac_z = (cptac_f - mu_c) / s_c
    print(f"  CPTAC: {cptac_z.shape[0]} samples × {cptac_z.shape[1]} genes")

    pca_cptac = PCA(n_components=10)
    pcs_cptac = pca_cptac.fit_transform(cptac_z.values)
    pve_cptac = pca_cptac.explained_variance_ratio_
    r2_cptac = tissue_r2(pcs_cptac, cptac_tissues)
    print(f"  PVE: {[f'{v*100:.1f}%' for v in pve_cptac[:5]]}")
    print(f"  Tissue R² per PC: {[f'{v:.4f}' for v in r2_cptac[:5]]}")
    print(f"  Total tissue R² (PC1-5): {sum(r2_cptac[:5] * pve_cptac[:5]) / sum(pve_cptac[:5]):.4f}")

    make_scatter(pcs_cptac, cptac_tissues, pve_cptac, "CPTAC Z-scored PCA (4 studies)",
                 OUTDIR / "cptac_pca_zscore.png")

    # ── Summary table ────────────────────────────────────────────────────
    summary = pd.DataFrame({
        "dataset": ["CCLE_raw", "CCLE_zscore", "CPTAC_zscore"],
        "n_samples": [len(ccle_f), len(ccle_z), len(cptac_z)],
        "n_genes": [ccle_f.shape[1], ccle_z.shape[1], cptac_z.shape[1]],
        "n_tissues": [len(set(tissues)), len(set(tissues)), len(set(cptac_tissues))],
        "pc1_tissue_r2": [r2_raw[0], r2_z[0], r2_cptac[0]],
        "pc2_tissue_r2": [r2_raw[1], r2_z[1], r2_cptac[1]],
        "pc1_pve": [pve_raw[0], pve_z[0], pve_cptac[0]],
        "weighted_tissue_r2_pc1to5": [
            sum(r2_raw[:5] * pve_raw[:5]) / sum(pve_raw[:5]),
            sum(r2_z[:5] * pve_z[:5]) / sum(pve_z[:5]),
            sum(r2_cptac[:5] * pve_cptac[:5]) / sum(pve_cptac[:5]),
        ],
    })
    summary.to_csv(OUTDIR / "tissue_structure_summary.csv", index=False)
    print(f"\n{'='*70}")
    print(summary.to_string(index=False))
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
