"""Task-specific sample subset construction.

Wraps existing R subset_strategies.R for complex mixture-balanced logic,
and implements simpler Python-native subsets where possible.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from harmonize.utils.r_bridge import run_r_script
from harmonize.utils.paths import ProjectPaths

logger = logging.getLogger(__name__)


def build_subtype_subset(
    paths: ProjectPaths,
    task_config: dict[str, Any],
) -> pd.DataFrame:
    """
    Build the breast subtype sample subset.

    Uses the R subset_strategies.R for CPTAC mixture-balanced selection
    and explicit CCLE cell line lists from config.
    """
    from harmonize.preprocessing.metadata import build_sample_meta
    from harmonize.preprocessing.loaders import load_cptac_studies, load_ccle

    cptac = load_cptac_studies(paths, {
        task_config["cptac"]["studies"][0]: {
            "tissue": "Breast",
            "path": paths.cptac_breast_matrix,
        }
    })
    ccle = load_ccle(paths)
    meta = build_sample_meta(cptac, ccle, task_config, None, repo_root=paths.root)

    logger.info("Subtype subset: %d samples (%d CPTAC, %d CCLE)",
                len(meta),
                (meta["domain"] == "CPTAC").sum(),
                (meta["domain"] == "CCLE").sum())
    return meta


def build_bvl_subset(
    paths: ProjectPaths,
    task_config: dict[str, Any],
) -> pd.DataFrame:
    """Build the breast vs lung sample subset."""
    from harmonize.preprocessing.metadata import build_sample_meta
    from harmonize.preprocessing.loaders import load_cptac_studies, load_ccle, load_ccle_sample_info

    studies = {}
    breast_study = task_config["cptac"].get("breast_study", "PDC000120")
    lung_study = task_config["cptac"].get("lung_study", "PDC000153")

    if paths.cptac_breast_matrix.exists():
        studies[breast_study] = {"tissue": "Breast", "path": paths.cptac_breast_matrix}
    if paths.cptac_lung_matrix.exists():
        studies[lung_study] = {"tissue": "Lung", "path": paths.cptac_lung_matrix}

    cptac = load_cptac_studies(paths, studies)
    ccle = load_ccle(paths)
    ccle_info = load_ccle_sample_info(paths)
    meta = build_sample_meta(cptac, ccle, task_config, ccle_info, repo_root=paths.root)

    logger.info("BvL subset: %d samples (%d CPTAC, %d CCLE)",
                len(meta),
                (meta["domain"] == "CPTAC").sum(),
                (meta["domain"] == "CCLE").sum())
    return meta
