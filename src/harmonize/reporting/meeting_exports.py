"""Collect key benchmark outputs into a meeting-ready export folder."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def export_meeting_materials(
    benchmark_root: Path,
    outdir: Path,
    config: dict[str, Any] | None = None,
) -> list[str]:
    """
    Gather the most important benchmark outputs into a meeting folder.

    Collects:
    - Polished marker profiles
    - Structure PCA/UMAP plots
    - Comparison summary tables
    - Meeting notes
    """
    outdir.mkdir(parents=True, exist_ok=True)
    exported = []

    # Polished profiles
    profiles_src = benchmark_root / "marker_profiles"
    for task in ["breast_subtype", "breast_vs_lung"]:
        polished = profiles_src / task / "polished"
        if polished.exists():
            dest = outdir / "profiles" / task
            dest.mkdir(parents=True, exist_ok=True)
            for f in polished.glob("*.png"):
                shutil.copy2(f, dest / f.name)
                exported.append(str(dest / f.name))

    # Structure plots
    struct_src = benchmark_root / "structure_metrics"
    if struct_src.exists():
        dest = outdir / "structure"
        dest.mkdir(parents=True, exist_ok=True)
        for f in struct_src.rglob("*.png"):
            shutil.copy2(f, dest / f.name)
            exported.append(str(dest / f.name))

    # Summary tables
    results_src = benchmark_root / "benchmark_results"
    if results_src.exists():
        dest = outdir / "tables"
        dest.mkdir(parents=True, exist_ok=True)
        for f in results_src.glob("*.csv"):
            shutil.copy2(f, dest / f.name)
            exported.append(str(dest / f.name))

    # Meeting notes
    meeting_src = benchmark_root / "meeting"
    if meeting_src.exists():
        for f in meeting_src.glob("*.md"):
            shutil.copy2(f, outdir / f.name)
            exported.append(str(outdir / f.name))

    # Comparison tables from diagnostics
    diag_comp = benchmark_root / "diagnostics" / "benchmark_comparison"
    if diag_comp.exists():
        dest = outdir / "tables"
        dest.mkdir(parents=True, exist_ok=True)
        for f in diag_comp.glob("*.csv"):
            shutil.copy2(f, dest / f.name)
            exported.append(str(dest / f.name))

    logger.info("Exported %d files to %s", len(exported), outdir)
    return exported
