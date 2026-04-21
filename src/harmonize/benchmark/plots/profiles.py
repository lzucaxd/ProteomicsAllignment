"""Wrapper for polished marker profile plots (calls existing R code)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from harmonize.utils.r_bridge import run_r_script
from harmonize.utils.paths import ProjectPaths

logger = logging.getLogger(__name__)


def generate_profile_plots(
    paths: ProjectPaths,
    config: dict[str, Any] | None = None,
) -> list[str]:
    """
    Generate polished marker profile plots by calling the existing R script.

    This wraps scripts/benchmark/run_polished_profile_plots.R which sources
    polished_profile_plots.R and generates PNG/PDF figures.
    """
    script = paths.r_benchmark_dir / "run_polished_profile_plots.R"
    if not script.exists():
        logger.warning("Profile plot script not found: %s", script)
        return []

    logger.info("Generating polished marker profile plots...")
    run_r_script(script, cwd=paths.root)

    # Collect generated files
    profile_dir = paths.resolve("reports/benchmark_master/marker_profiles")
    generated = []
    for task_dir in ["breast_subtype", "breast_vs_lung"]:
        polished = profile_dir / task_dir / "polished"
        if polished.exists():
            for f in polished.glob("*.png"):
                generated.append(str(f))

    logger.info("Generated %d profile plot files", len(generated))
    return generated
