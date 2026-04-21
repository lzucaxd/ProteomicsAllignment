#!/usr/bin/env python3
"""
Full Celligner alignment: all CPTAC studies × CCLE.

Pipeline:
  1. Union of CPTAC gene matrices (outer join on genes)
  2. Intersect with CCLE genes
  3. 70% prevalence filter per domain
  4. Remove near-constant genes (SD < 0.01 in either domain)
  5. Within-domain gene-median imputation
  6. Z-score standardization (per gene, per domain)
  7. Celligner: fit on CCLE, transform CPTAC
  8. PCA + UMAP before and after, colored by primary site
  9. Re-run benchmark
"""
import os, sys, time, json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "models" / "celligner-master"))

CPTAC_STUDIES = {
    "PDC000120": ("Breast", REPO / "data/results/PDC000120/gene_matrix.csv"),
    "PDC000127": ("Ovarian", REPO / "data/results/PDC000127/gene_matrix.csv"),
    "PDC000153": ("Lung", REPO / "data/results/PDC000153/gene_matrix.csv"),
    "PDC000204": ("Uterine", REPO / "data/results/PDC000204/gene_matrix.csv"),
}
CCLE_PATH = REPO / "data/results/CCLE_corrected/gene_matrix.csv"
CCLE_INFO = REPO / "data/ccle_peptide/sample_info_ccle.csv"
OUTDIR = REPO / "reports/benchmark_master/celligner_all"

MIN_PREVALENCE = 0.70
MIN_SD = 0.01


def load_gm(path):
    df = pd.read_csv(path, index_col=0)
    if "UniProtID" in df.columns:
        df = df.drop(columns=["UniProtID"])
    return df.T.astype(float)  # samples × genes


def get_ccle_tissue_map():
    info = pd.read_csv(CCLE_INFO)
    return dict(zip(info["Cell Line"], info["Tissue of Origin"]))


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    print("=" * 70)
    print("  CELLIGNER ALL-DATA ALIGNMENT")
    print("=" * 70)

    # ── Step 1: Load and union CPTAC ─────────────────────────────────────
    print("\n[1/8] Loading CPTAC studies (outer join on genes)...")
    cptac_parts = {}
    cptac_study_labels = {}
    cptac_cancer_labels = {}

    for study_id, (cancer, path) in CPTAC_STUDIES.items():
        if not path.exists():
            print(f"  SKIP {study_id}: file not found")
            continue
        df = load_gm(path)
        print(f"  {study_id} ({cancer}): {df.shape[0]} samples × {df.shape[1]} genes")
        cptac_parts[study_id] = df
        for s in df.index:
            cptac_study_labels[s] = study_id
            cptac_cancer_labels[s] = cancer

    cptac_union = pd.concat(cptac_parts.values(), axis=0, join="outer")
    print(f"  CPTAC union: {cptac_union.shape[0]} samples × {cptac_union.shape[1]} genes")

    # ── Step 2: Load CCLE and intersect ──────────────────────────────────
    print("\n[2/8] Loading CCLE and intersecting genes...")
    ccle = load_gm(CCLE_PATH)
    print(f"  CCLE: {ccle.shape[0]} samples × {ccle.shape[1]} genes")

    shared_genes = sorted(set(cptac_union.columns) & set(ccle.columns))
    print(f"  Shared genes (CPTAC union ∩ CCLE): {len(shared_genes)}")
    cptac_s = cptac_union[shared_genes]
    ccle_s = ccle[shared_genes]

    # ── Step 3: 70% prevalence filter per domain ─────────────────────────
    print(f"\n[3/8] Prevalence filter ({MIN_PREVALENCE*100:.0f}% in each domain)...")
    cptac_prev = cptac_s.notna().mean(axis=0)
    ccle_prev = ccle_s.notna().mean(axis=0)
    prev_mask = (cptac_prev >= MIN_PREVALENCE) & (ccle_prev >= MIN_PREVALENCE)
    genes_prev = [g for g in shared_genes if prev_mask[g]]
    print(f"  Genes passing prevalence: {len(genes_prev)} (dropped {len(shared_genes) - len(genes_prev)})")
    cptac_s = cptac_s[genes_prev]
    ccle_s = ccle_s[genes_prev]

    # ── Step 4: Remove near-constant genes ───────────────────────────────
    print(f"\n[4/8] Removing near-constant genes (SD < {MIN_SD})...")
    cptac_sd = cptac_s.std(axis=0, skipna=True)
    ccle_sd = ccle_s.std(axis=0, skipna=True)
    sd_mask = (cptac_sd >= MIN_SD) & (ccle_sd >= MIN_SD)
    genes_sd = [g for g in genes_prev if sd_mask[g]]
    n_dropped = len(genes_prev) - len(genes_sd)
    print(f"  Genes passing SD filter: {len(genes_sd)} (dropped {n_dropped})")
    cptac_s = cptac_s[genes_sd]
    ccle_s = ccle_s[genes_sd]

    # ── Step 5: Within-domain gene-median imputation ─────────────────────
    print("\n[5/8] Within-domain gene-median imputation...")
    na_cptac = cptac_s.isna().sum().sum()
    na_ccle = ccle_s.isna().sum().sum()
    cptac_s = cptac_s.fillna(cptac_s.median(axis=0))
    ccle_s = ccle_s.fillna(ccle_s.median(axis=0))
    print(f"  CPTAC NAs imputed: {na_cptac}")
    print(f"  CCLE NAs imputed:  {na_ccle}")

    # ── Step 6: Z-score standardization (per gene, per domain) ───────────
    print("\n[6/8] Z-score standardization (per gene, per domain)...")
    cptac_mean = cptac_s.mean(axis=0)
    cptac_std = cptac_s.std(axis=0)
    cptac_std[cptac_std == 0] = 1.0
    cptac_z = (cptac_s - cptac_mean) / cptac_std

    ccle_mean = ccle_s.mean(axis=0)
    ccle_std = ccle_s.std(axis=0)
    ccle_std[ccle_std == 0] = 1.0
    ccle_z = (ccle_s - ccle_mean) / ccle_std

    print(f"  CPTAC: {cptac_z.shape}")
    print(f"  CCLE:  {ccle_z.shape}")

    n_genes = len(genes_sd)
    n_cptac = len(cptac_z)
    n_ccle = len(ccle_z)

    # ── Build metadata ───────────────────────────────────────────────────
    tissue_map = get_ccle_tissue_map()
    tissue_norm_map = {k.replace("-", "").replace("_", "").upper(): v for k, v in tissue_map.items()}

    def get_site(sample_id):
        if sample_id in cptac_cancer_labels:
            return cptac_cancer_labels[sample_id]
        if sample_id in tissue_map:
            return tissue_map[sample_id]
        s_norm = sample_id.replace("-", "").replace("_", "").upper()
        for k, v in tissue_norm_map.items():
            if k == s_norm:
                return v
        return "Unknown"

    def get_domain(s):
        return "CPTAC" if s in cptac_study_labels else "CCLE"

    # ── Pre-alignment PCA + UMAP ─────────────────────────────────────────
    print("\n[7/8] Running Celligner...")

    combined_pre = pd.concat([cptac_z, ccle_z])
    domains_all = [get_domain(s) for s in combined_pre.index]
    sites_all = [get_site(s) for s in combined_pre.index]

    # PCA pre
    pca_pre = PCA(n_components=2)
    pcs_pre = pca_pre.fit_transform(combined_pre.values)
    pve_pre = pca_pre.explained_variance_ratio_
    d_num = np.array([0 if d == "CPTAC" else 1 for d in domains_all]).reshape(-1, 1)
    r2_pre = LinearRegression().fit(d_num, pcs_pre[:, 0]).score(d_num, pcs_pre[:, 0])
    print(f"  Pre-alignment PC1 domain R² = {r2_pre:.4f}, PVE = {pve_pre[0]*100:.1f}%")

    # Run Celligner
    from celligner import Celligner
    model = Celligner()
    print(f"  Fitting on CCLE ({ccle_z.shape})...")
    model.fit(ccle_z)
    print(f"  Transforming CPTAC ({cptac_z.shape})...")
    model.transform(cptac_z)
    elapsed_cell = time.time() - t_start
    combined_post = model.combined_output
    print(f"  Celligner done in {elapsed_cell:.1f}s — output: {combined_post.shape}")

    # MNN pairs
    mnn_count = "unknown"
    if hasattr(model, "mnn_pairs"):
        mnn_count = len(model.mnn_pairs) if model.mnn_pairs is not None else "none"

    # Post-alignment PCA
    domains_post = [get_domain(s) for s in combined_post.index]
    sites_post = [get_site(s) for s in combined_post.index]

    pca_post = PCA(n_components=2)
    pcs_post = pca_post.fit_transform(combined_post.values)
    pve_post = pca_post.explained_variance_ratio_
    d_num2 = np.array([0 if d == "CPTAC" else 1 for d in domains_post]).reshape(-1, 1)
    r2_post = LinearRegression().fit(d_num2, pcs_post[:, 0]).score(d_num2, pcs_post[:, 0])
    print(f"  Post-alignment PC1 domain R² = {r2_post:.4f}, PVE = {pve_post[0]*100:.1f}%")

    # ── Step 8: Plots ────────────────────────────────────────────────────
    print("\n[8/8] Generating PCA + UMAP plots...")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # Collapse tissue labels
    TISSUE_REMAP = {
        "Breast": "Breast", "Lung": "Lung", "Ovarian": "Ovary", "Uterine": "Uterine",
        "Ovary": "Ovary", "Large Intestine": "Colorectal",
        "Haematopoietic and Lymphoid Tissue": "Blood/Lymphoid",
        "Acute Myeloid Leukemia": "Blood/Lymphoid", "Lymphoma": "Blood/Lymphoid",
        "Skin": "Skin", "Central Nervous System": "CNS",
        "Pancreas": "Pancreas", "Stomach": "Upper GI", "Oesophagus": "Upper GI",
        "Endometrium": "Endometrium", "Liver": "Liver", "Kidney": "Kidney",
        "Urinary Tract": "Urinary Tract", "Upper Aerodigestive Tract": "Head & Neck",
        "Prostate": "Prostate", "Soft Tissue": "Other", "Bone": "Other",
        "Thyroid": "Other", "Pleura": "Other", "Autonomic Ganglia": "Other",
        "Biliary Tract": "Other",
    }
    SITE_COLORS = {
        "Breast": "#E91E63", "Lung": "#4CAF50", "Ovary": "#00BCD4",
        "Uterine": "#F06292", "Colorectal": "#FF9800",
        "Blood/Lymphoid": "#9C27B0", "Skin": "#795548", "CNS": "#3F51B5",
        "Pancreas": "#CDDC39", "Upper GI": "#FF5722", "Endometrium": "#AB47BC",
        "Liver": "#8D6E63", "Kidney": "#5C6BC0", "Urinary Tract": "#26A69A",
        "Head & Neck": "#78909C", "Prostate": "#AED581", "Other": "#BDBDBD",
        "Unknown": "#E0E0E0",
    }
    DOMAIN_MARKERS = {"CPTAC": "^", "CCLE": "o"}

    def remap_site(s):
        return TISSUE_REMAP.get(s, "Other")

    def scatter_plot(coords, sites, domains, pve, r2, title, outpath):
        fig, ax = plt.subplots(figsize=(14, 10))
        sites_r = [remap_site(s) for s in sites]
        for site in SITE_COLORS:
            for dom in ("CCLE", "CPTAC"):
                idx = [i for i, (s, d) in enumerate(zip(sites_r, domains)) if s == site and d == dom]
                if not idx: continue
                ax.scatter(coords[idx, 0], coords[idx, 1],
                           c=SITE_COLORS[site], marker=DOMAIN_MARKERS[dom],
                           alpha=0.7, s=30, edgecolors="k" if dom == "CPTAC" else "none",
                           linewidths=0.3)
        from collections import Counter
        ct = Counter(sites_r)
        present = sorted(ct.keys(), key=lambda s: -ct[s])
        handles = [Line2D([0],[0], marker='o', color='w',
                   markerfacecolor=SITE_COLORS.get(s,"gray"), markersize=7,
                   label=f"{s} (n={ct[s]})") for s in present]
        handles += [Line2D([0],[0], marker='^', color='w', markerfacecolor='gray',
                    markeredgecolor='k', markersize=7, label='CPTAC'),
                    Line2D([0],[0], marker='o', color='w', markerfacecolor='gray',
                    markersize=7, label='CCLE')]
        leg = ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left",
                        fontsize=7, title="Primary Site / Domain", title_fontsize=8)
        ax.set_xlabel(f"Dim 1 ({pve[0]*100:.1f}%)" if pve is not None else "Dim 1")
        ax.set_ylabel(f"Dim 2 ({pve[1]*100:.1f}%)" if pve is not None else "Dim 2")
        ax.set_title(f"{title}  (domain R²={r2:.4f})")
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        plt.close()
        print(f"  Saved: {outpath.name}")

    # PCA plots
    scatter_plot(pcs_pre, sites_all, domains_all, pve_pre, r2_pre,
                 "Pre-Alignment PCA", OUTDIR / "pca_pre.png")
    scatter_plot(pcs_post, sites_post, domains_post, pve_post, r2_post,
                 "Post-Alignment PCA", OUTDIR / "pca_post.png")

    # UMAP
    print("  Computing UMAP (pre)...")
    import umap
    um_pre = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=42, n_jobs=1)
    umap_pre = um_pre.fit_transform(combined_pre.values)
    scatter_plot(umap_pre, sites_all, domains_all, None, r2_pre,
                 "Pre-Alignment UMAP", OUTDIR / "umap_pre.png")

    print("  Computing UMAP (post)...")
    um_post = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=42, n_jobs=1)
    umap_post = um_post.fit_transform(combined_post.values)
    scatter_plot(umap_post, sites_post, domains_post, None, r2_post,
                 "Post-Alignment UMAP", OUTDIR / "umap_post.png")

    # Save aligned matrix (genes × samples for benchmark)
    print("\n  Saving outputs...")
    combined_post.to_csv(OUTDIR / "celligner_aligned_matrix.csv")

    # Save metadata
    meta_rows = []
    for s in combined_post.index:
        meta_rows.append({
            "sample_id": s,
            "domain": get_domain(s),
            "primary_site": get_site(s),
            "primary_site_collapsed": remap_site(get_site(s)),
            "study_id": cptac_study_labels.get(s, "CCLE"),
        })
    pd.DataFrame(meta_rows).to_csv(OUTDIR / "sample_metadata.csv", index=False)

    # Save PCA/UMAP coords
    pd.DataFrame({"PC1": pcs_post[:, 0], "PC2": pcs_post[:, 1],
                   "UMAP1": umap_post[:, 0], "UMAP2": umap_post[:, 1],
                   "domain": domains_post, "site": sites_post,
                   "sample": combined_post.index}).to_csv(OUTDIR / "coords_post.csv", index=False)

    # Save DE genes
    if hasattr(model, "de_genes") and model.de_genes is not None:
        with open(OUTDIR / "de_genes.txt", "w") as f:
            f.write("\n".join(model.de_genes))

    summary = {
        "cptac_studies": list(CPTAC_STUDIES.keys()),
        "cptac_samples": int(n_cptac),
        "ccle_samples": int(n_ccle),
        "total_samples": int(n_cptac + n_ccle),
        "genes_shared": int(len(shared_genes)),
        "genes_after_prevalence": int(len(genes_prev)),
        "genes_after_sd": int(len(genes_sd)),
        "prevalence_threshold": MIN_PREVALENCE,
        "sd_threshold": MIN_SD,
        "pc1_domain_r2_pre": float(round(r2_pre, 6)),
        "pc1_domain_r2_post": float(round(r2_post, 6)),
        "runtime_seconds": float(round(time.time() - t_start, 1)),
    }
    with open(OUTDIR / "run_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*70}")
    print(f"  ALL-DATA CELLIGNER COMPLETE")
    print(f"  Samples: {n_cptac} CPTAC + {n_ccle} CCLE = {n_cptac+n_ccle}")
    print(f"  Genes: {len(shared_genes)} shared → {len(genes_prev)} (prevalence) → {len(genes_sd)} (SD)")
    print(f"  Domain R²: {r2_pre:.4f} → {r2_post:.4f}")
    print(f"  Runtime: {time.time()-t_start:.1f}s")
    print(f"  Outputs: {OUTDIR}/")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
