"""Load CPTAC and CCLE gene matrices from processed results."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from harmonize.utils.io import load_gene_matrix
from harmonize.utils.paths import ProjectPaths

logger = logging.getLogger(__name__)


def load_cptac_studies(
    paths: ProjectPaths, studies: dict[str, dict[str, Any]] | None = None
) -> dict[str, pd.DataFrame]:
    """
    Load CPTAC gene matrices for all configured studies.

    Returns dict mapping study_id -> DataFrame (genes x samples).
    """
    if studies is None:
        studies = {
            "PDC000120": {"tissue": "Breast", "path": paths.cptac_breast_matrix},
            "PDC000153": {"tissue": "Lung", "path": paths.cptac_lung_matrix},
        }

    result = {}
    for study_id, info in studies.items():
        mat_path = info.get("path") or info.get("gene_matrix")
        if mat_path is None:
            continue
        mat_path = paths.resolve(str(mat_path)) if not Path(str(mat_path)).is_absolute() else Path(mat_path)
        if not mat_path.exists():
            logger.warning("CPTAC %s matrix not found: %s", study_id, mat_path)
            continue
        df = load_gene_matrix(mat_path)
        logger.info("Loaded CPTAC %s (%s): %d genes x %d samples", study_id, info.get("tissue", "?"), len(df), len(df.columns))
        result[study_id] = df

    return result


def load_ccle(paths: ProjectPaths) -> pd.DataFrame:
    """Load the corrected CCLE gene matrix."""
    df = load_gene_matrix(paths.ccle_matrix)
    logger.info("Loaded CCLE: %d genes x %d samples", len(df), len(df.columns))
    return df


def load_ccle_sample_info(paths: ProjectPaths) -> pd.DataFrame:
    """Load CCLE sample info with cell line and tissue annotations."""
    df = pd.read_csv(paths.ccle_sample_info, on_bad_lines="skip")
    col_renames = {}
    for c in df.columns:
        if c.replace(".", " ") == "Cell Line":
            col_renames[c] = "Cell Line"
        elif c.replace(".", " ") == "Tissue of Origin":
            col_renames[c] = "Tissue of Origin"
    df = df.rename(columns=col_renames)
    df = df[df["Cell Line"].notna() & (df["Cell Line"].str.len() > 0)]
    return df
