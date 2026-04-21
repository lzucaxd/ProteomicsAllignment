"""Diagnostic plot generation: QQ plots, distribution checks, effect summaries."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from harmonize.utils.r_bridge import run_r_script
from harmonize.utils.paths import ProjectPaths

logger = logging.getLogger(__name__)


def generate_diagnostics(
    paths: ProjectPaths,
    config: dict[str, Any] | None = None,
) -> list[str]:
    """Run existing R diagnostics script."""
    script = paths.r_benchmark_dir / "diagnostics.R"
    if not script.exists():
        logger.warning("Diagnostics R script not found: %s", script)
        return []

    run_r_script(script, cwd=paths.root)
    diag_dir = paths.resolve("reports/benchmark_master/diagnostics")
    return [str(f) for f in diag_dir.rglob("*.png")]


def plot_fc_scatter(
    da_cptac: pd.DataFrame,
    da_ccle: pd.DataFrame,
    method_name: str,
    task_name: str,
    outdir: str | Path,
    gene_col: str = "gene",
    fc_col: str = "logFC",
    markers: list[str] | None = None,
) -> str | None:
    """Fold-change scatter: CPTAC logFC vs CCLE logFC."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    merged = pd.merge(
        da_cptac[[gene_col, fc_col]].rename(columns={fc_col: "fc_cptac"}),
        da_ccle[[gene_col, fc_col]].rename(columns={fc_col: "fc_ccle"}),
        on=gene_col,
    )
    if len(merged) < 5:
        return None

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(merged["fc_cptac"], merged["fc_ccle"], alpha=0.15, s=8, c="gray")

    if markers:
        mk = merged[merged[gene_col].isin(markers)]
        ax.scatter(mk["fc_cptac"], mk["fc_ccle"], c="red", s=30, zorder=5, edgecolors="black", linewidths=0.5)
        for _, row in mk.iterrows():
            ax.annotate(row[gene_col], (row["fc_cptac"], row["fc_ccle"]),
                        fontsize=7, ha="left", va="bottom")

    lim = max(abs(merged["fc_cptac"]).max(), abs(merged["fc_ccle"]).max()) * 1.1
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axhline(0, c="gray", lw=0.5, ls="--")
    ax.axvline(0, c="gray", lw=0.5, ls="--")
    ax.plot([-lim, lim], [-lim, lim], c="blue", lw=0.8, ls=":", alpha=0.5)
    ax.set_xlabel("CPTAC logFC")
    ax.set_ylabel("CCLE logFC")
    ax.set_title(f"{method_name} — {task_name}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    corr = merged[["fc_cptac", "fc_ccle"]].corr().iloc[0, 1]
    same = (np.sign(merged["fc_cptac"]) == np.sign(merged["fc_ccle"])).mean()
    ax.text(0.05, 0.95, f"r = {corr:.3f}\nsame dir = {same:.1%}",
            transform=ax.transAxes, va="top", fontsize=9)

    p = outdir / f"fc_scatter_{method_name}_{task_name}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(p)
