"""Structure visualization: PCA and UMAP colored by domain/condition."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

logger = logging.getLogger(__name__)

DOMAIN_COLORS = {"CPTAC": "#1565C0", "CCLE": "#E53935"}


def plot_pca_structure(
    matrix: pd.DataFrame,
    sample_meta: pd.DataFrame,
    method_name: str,
    task_name: str,
    outdir: str | Path,
    n_pcs: int = 2,
) -> list[str]:
    """Generate PCA scatter plots colored by domain and condition."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = sample_meta.set_index("sample_id") if "sample_id" in sample_meta.columns else sample_meta
    valid = [c for c in matrix.columns if c in meta.index]
    if len(valid) < 10:
        logger.warning("Too few samples (%d) for PCA plot", len(valid))
        return []

    X = matrix[valid].T.fillna(0).values
    meta = meta.loc[valid]

    pca = PCA(n_components=min(n_pcs, min(X.shape) - 1))
    Z = pca.fit_transform(X)
    var_exp = pca.explained_variance_ratio_

    paths_out = []

    # Plot by domain
    if "domain" in meta.columns:
        fig, ax = plt.subplots(figsize=(7, 5.5))
        for domain in meta["domain"].unique():
            mask = meta["domain"].values == domain
            color = DOMAIN_COLORS.get(domain, "#888888")
            ax.scatter(Z[mask, 0], Z[mask, 1], c=color, label=domain,
                       alpha=0.5, s=15, edgecolors="none")
        ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)")
        ax.set_title(f"{method_name} — {task_name} (by domain)")
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        p = outdir / f"pca_domain_{method_name}_{task_name}.png"
        fig.savefig(p, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths_out.append(str(p))

    # Plot by condition
    if "condition" in meta.columns:
        fig, ax = plt.subplots(figsize=(7, 5.5))
        conditions = meta["condition"].unique()
        cmap = plt.cm.Set2(np.linspace(0, 1, max(len(conditions), 3)))
        for i, cond in enumerate(conditions):
            mask = meta["condition"].values == cond
            ax.scatter(Z[mask, 0], Z[mask, 1], c=[cmap[i]], label=cond,
                       alpha=0.5, s=15, edgecolors="none")
        ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)")
        ax.set_title(f"{method_name} — {task_name} (by condition)")
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        p = outdir / f"pca_condition_{method_name}_{task_name}.png"
        fig.savefig(p, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths_out.append(str(p))

    return paths_out


def plot_umap_structure(
    matrix: pd.DataFrame,
    sample_meta: pd.DataFrame,
    method_name: str,
    task_name: str,
    outdir: str | Path,
    n_neighbors: int = 15,
    min_dist: float = 0.3,
) -> list[str]:
    """Generate UMAP scatter plots (requires umap-learn)."""
    try:
        import umap
    except ImportError:
        logger.warning("umap-learn not installed, skipping UMAP plots")
        return []

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = sample_meta.set_index("sample_id") if "sample_id" in sample_meta.columns else sample_meta
    valid = [c for c in matrix.columns if c in meta.index]
    if len(valid) < 15:
        return []

    X = matrix[valid].T.fillna(0).values
    meta = meta.loc[valid]

    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, random_state=42)
    Z = reducer.fit_transform(X)

    paths_out = []

    if "domain" in meta.columns:
        fig, ax = plt.subplots(figsize=(7, 5.5))
        for domain in meta["domain"].unique():
            mask = meta["domain"].values == domain
            color = DOMAIN_COLORS.get(domain, "#888888")
            ax.scatter(Z[mask, 0], Z[mask, 1], c=color, label=domain,
                       alpha=0.5, s=15, edgecolors="none")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"{method_name} — {task_name} (UMAP, by domain)")
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        p = outdir / f"umap_domain_{method_name}_{task_name}.png"
        fig.savefig(p, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths_out.append(str(p))

    return paths_out
