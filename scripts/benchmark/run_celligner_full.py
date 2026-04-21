#!/usr/bin/env python3
"""
Full Celligner run: CPTAC breast (PDC000120) + CCLE (all cell lines).

Usage:
  R_HOME=/Library/Frameworks/R.framework/Resources \
  R_PROFILE_USER=/dev/null \
  python scripts/benchmark/run_celligner_full.py
"""
import os, sys, time, json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "models" / "celligner-master"))

CPTAC_BREAST = REPO / "data" / "results" / "PDC000120" / "gene_matrix.csv"
CCLE_MATRIX  = REPO / "data" / "results" / "CCLE_corrected" / "gene_matrix.csv"
OUTDIR       = REPO / "reports" / "benchmark_master" / "celligner_full"

MIN_OBS_FRAC = 0.5

SUBTYPE_BASAL   = ["HCC70", "HCC1806", "HCC1143", "MDA-MB-468"]
SUBTYPE_LUMINAL = ["CAMA-1", "MCF7", "T-47D", "ZR-75-1"]

MARKERS_SUBTYPE = {
    "FOXA1": +1, "GATA3": +1, "ESR1": +1, "PGR": +1,
    "KRT5": -1, "KRT14": -1, "KRT17": -1, "EGFR": -1,
    "ERBB2": +1, "CDH1": +1,
}

CCLE_SAMPLE_INFO = REPO / "data" / "ccle_peptide" / "sample_info_ccle.csv"


def load_matrix(path):
    df = pd.read_csv(path, index_col=0)
    if "UniProtID" in df.columns:
        df = df.drop(columns=["UniProtID"])
    return df.T.astype(float)


def pc1_domain_r2(mat, labels):
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(mat)
    d = np.array([0 if l == "CPTAC" else 1 for l in labels]).reshape(-1, 1)
    r2 = LinearRegression().fit(d, pcs[:, 0]).score(d, pcs[:, 0])
    return r2, pca.explained_variance_ratio_[:2], pcs


def match_subtype_lines(index, line_names):
    out = []
    norm = {c.replace("-", "").replace("_", "").upper(): c for c in index}
    for ln in line_names:
        ln_norm = ln.replace("-", "").replace("_", "").upper()
        for k, v in norm.items():
            if ln_norm in k:
                out.append(v)
                break
    return out


def make_pca_plot(pcs, labels, color_map, title, outpath, shape_col=None, shapes=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 8))
    unique_labels = sorted(set(labels))

    if shape_col is not None:
        unique_shapes = sorted(set(shape_col))
        shape_markers = {"CPTAC": "o", "CCLE": "^"}
        for lab in unique_labels:
            for sh in unique_shapes:
                mask = [(l == lab and s == sh) for l, s in zip(labels, shape_col)]
                if not any(mask):
                    continue
                idx = [i for i, m in enumerate(mask) if m]
                ax.scatter(pcs[idx, 0], pcs[idx, 1],
                           c=color_map.get(lab, "gray"),
                           marker=shape_markers.get(sh, "o"),
                           label=f"{lab} ({sh})", alpha=0.7, s=40, edgecolors="k", linewidths=0.3)
    else:
        for lab in unique_labels:
            idx = [i for i, l in enumerate(labels) if l == lab]
            ax.scatter(pcs[idx, 0], pcs[idx, 1],
                       c=color_map.get(lab, "gray"),
                       label=lab, alpha=0.7, s=40, edgecolors="k", linewidths=0.3)

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  Plot saved: {outpath}")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  CELLIGNER FULL RUN — CPTAC breast + CCLE")
    print("=" * 60)

    # Load
    print("\n[1/6] Loading matrices...")
    cptac = load_matrix(CPTAC_BREAST)
    ccle  = load_matrix(CCLE_MATRIX)
    print(f"  CPTAC: {cptac.shape}")
    print(f"  CCLE:  {ccle.shape}")

    # Intersect and clean
    shared = sorted(set(cptac.columns) & set(ccle.columns))
    cptac_s = cptac[shared]
    ccle_s  = ccle[shared]
    obs = (cptac_s.notna().mean(0) >= MIN_OBS_FRAC) & (ccle_s.notna().mean(0) >= MIN_OBS_FRAC)
    genes = [g for g in shared if obs[g]]
    cptac_s = cptac_s[genes].fillna(cptac_s[genes].median(0))
    ccle_s  = ccle_s[genes].fillna(ccle_s[genes].median(0))
    print(f"  Shared genes after filter: {len(genes)}")

    # Pre-alignment
    print("\n[2/6] Pre-alignment metrics...")
    combined_pre = pd.concat([cptac_s, ccle_s])
    domains_pre = ["CPTAC"] * len(cptac_s) + ["CCLE"] * len(ccle_s)
    r2_pre, pve_pre, pcs_pre = pc1_domain_r2(combined_pre.values, domains_pre)
    print(f"  PC1 domain R²: {r2_pre:.4f}")
    print(f"  PC1 PVE: {pve_pre[0]*100:.1f}%")

    # Run Celligner
    print("\n[3/6] Running Celligner...")
    t0 = time.time()
    from celligner import Celligner
    model = Celligner()
    print(f"  Fitting on CCLE ({ccle_s.shape})...")
    model.fit(ccle_s)
    print(f"  Transforming CPTAC ({cptac_s.shape})...")
    model.transform(cptac_s)
    elapsed = time.time() - t0
    combined_post = model.combined_output
    print(f"  Done in {elapsed:.1f}s — output shape: {combined_post.shape}")

    # Post-alignment
    print("\n[4/6] Post-alignment metrics...")
    cptac_ids = set(cptac_s.index)
    domains_post = ["CPTAC" if s in cptac_ids else "CCLE" for s in combined_post.index]
    r2_post, pve_post, pcs_post = pc1_domain_r2(combined_post.values, domains_post)
    print(f"  PC1 domain R²: {r2_post:.4f}")
    print(f"  PC1 PVE: {pve_post[0]*100:.1f}%")
    print(f"  Domain R²: {r2_pre:.4f} → {r2_post:.4f} (Δ = {r2_post - r2_pre:+.4f})")

    # Marker check
    print("\n[5/6] Marker direction check...")
    ccle_aligned = combined_post.loc[[s for s in combined_post.index if s not in cptac_ids]]
    basal_ids = match_subtype_lines(ccle_aligned.index, SUBTYPE_BASAL)
    luminal_ids = match_subtype_lines(ccle_aligned.index, SUBTYPE_LUMINAL)
    print(f"  Basal lines found: {len(basal_ids)} — {basal_ids}")
    print(f"  Luminal lines found: {len(luminal_ids)} — {luminal_ids}")

    marker_results = []
    for gene, exp_sign in MARKERS_SUBTYPE.items():
        if gene not in combined_post.columns:
            marker_results.append({"gene": gene, "expected": exp_sign, "FC": np.nan, "status": "missing"})
            continue
        if not basal_ids or not luminal_ids:
            marker_results.append({"gene": gene, "expected": exp_sign, "FC": np.nan, "status": "no_samples"})
            continue
        lum = ccle_aligned.loc[luminal_ids, gene].mean()
        bas = ccle_aligned.loc[basal_ids, gene].mean()
        fc = lum - bas
        correct = "OK" if np.sign(fc) == exp_sign else "WRONG"
        marker_results.append({"gene": gene, "expected": exp_sign, "FC": fc, "status": correct})
        print(f"  {gene}: FC={fc:+.3f} expected={'+' if exp_sign>0 else '-'} → {correct}")

    n_ok = sum(1 for m in marker_results if m["status"] == "OK")
    n_chk = sum(1 for m in marker_results if m["status"] in ("OK", "WRONG"))
    print(f"  Direction accuracy: {n_ok}/{n_chk}")

    # Plots
    print("\n[6/6] Generating plots...")
    domain_colors = {"CPTAC": "#2196F3", "CCLE": "#FF5722"}

    make_pca_plot(pcs_pre, domains_pre, domain_colors,
                  f"Pre-Alignment PCA (PC1 domain R²={r2_pre:.3f})",
                  OUTDIR / "pca_pre_alignment.png")
    make_pca_plot(pcs_post, domains_post, domain_colors,
                  f"Post-Alignment PCA (PC1 domain R²={r2_post:.4f})",
                  OUTDIR / "pca_post_alignment.png")

    # Tissue-type PCA coloring (post-alignment)
    if CCLE_SAMPLE_INFO.exists():
        info = pd.read_csv(CCLE_SAMPLE_INFO)
        tissue_map = dict(zip(info["Cell Line"], info["Tissue of Origin"]))
        tissue_labels = []
        for s, d in zip(combined_post.index, domains_post):
            if d == "CPTAC":
                tissue_labels.append("CPTAC_Breast")
            else:
                t = tissue_map.get(s, "Other")
                tissue_labels.append(f"CCLE_{t}" if t in ("Breast", "Lung") else "CCLE_Other")
        tissue_colors = {
            "CPTAC_Breast": "#2196F3", "CCLE_Breast": "#E91E63",
            "CCLE_Lung": "#4CAF50", "CCLE_Other": "#BDBDBD"
        }
        make_pca_plot(pcs_post, tissue_labels, tissue_colors,
                      "Post-Alignment PCA by Tissue",
                      OUTDIR / "pca_post_tissue.png")

    # Save outputs
    combined_post.to_csv(OUTDIR / "celligner_aligned_matrix.csv")
    pd.DataFrame(marker_results).to_csv(OUTDIR / "marker_check.csv", index=False)

    pd.DataFrame({
        "PC1": pcs_post[:, 0], "PC2": pcs_post[:, 1],
        "domain": domains_post, "sample": combined_post.index
    }).to_csv(OUTDIR / "pca_post_data.csv", index=False)

    # Save DE genes
    if hasattr(model, "de_genes") and model.de_genes is not None:
        with open(OUTDIR / "celligner_de_genes.txt", "w") as f:
            f.write("\n".join(model.de_genes))
        print(f"  DE genes saved: {len(model.de_genes)}")

    summary = {
        "cptac_samples": int(len(cptac_s)),
        "ccle_samples": int(len(ccle_s)),
        "genes": int(len(genes)),
        "pc1_domain_r2_pre": float(round(r2_pre, 6)),
        "pc1_domain_r2_post": float(round(r2_post, 6)),
        "marker_direction_accuracy": f"{n_ok}/{n_chk}",
        "runtime_seconds": float(round(elapsed, 1)),
    }
    with open(OUTDIR / "run_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"  FULL CELLIGNER RUN COMPLETE")
    print(f"  Domain R²: {r2_pre:.4f} → {r2_post:.4f}")
    print(f"  Markers: {n_ok}/{n_chk} correct")
    print(f"  Runtime: {elapsed:.1f}s")
    print(f"  Outputs: {OUTDIR}/")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
