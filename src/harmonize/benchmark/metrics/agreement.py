"""Cross-domain fold-change agreement metrics."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_fc_agreement(
    da_cptac: pd.DataFrame,
    da_ccle: pd.DataFrame,
    fc_col: str = "logFC",
    gene_col: str = "gene",
) -> pd.DataFrame:
    """
    Compute fold-change agreement between CPTAC and CCLE DA results.

    Returns a merged DataFrame with columns:
    gene, logFC_cptac, logFC_ccle, same_direction, fc_diff
    """
    merged = pd.merge(
        da_cptac[[gene_col, fc_col]].rename(columns={fc_col: "logFC_cptac"}),
        da_ccle[[gene_col, fc_col]].rename(columns={fc_col: "logFC_ccle"}),
        on=gene_col,
        how="inner",
    )

    merged["same_direction"] = (
        np.sign(merged["logFC_cptac"]) == np.sign(merged["logFC_ccle"])
    )
    merged["fc_diff"] = merged["logFC_cptac"] - merged["logFC_ccle"]
    merged["fc_ratio"] = merged["logFC_cptac"] / merged["logFC_ccle"].replace(0, np.nan)

    return merged


def summarize_agreement(agreement: pd.DataFrame) -> dict:
    """Compute summary statistics from agreement DataFrame."""
    n = len(agreement)
    if n == 0:
        return {"n_genes": 0, "same_direction_frac": float("nan")}

    same_dir = agreement["same_direction"].sum()
    fc_corr = agreement[["logFC_cptac", "logFC_ccle"]].corr().iloc[0, 1]

    return {
        "n_genes": n,
        "same_direction_n": int(same_dir),
        "same_direction_frac": same_dir / n,
        "fc_correlation": fc_corr,
        "median_fc_diff": agreement["fc_diff"].median(),
        "mad_fc_diff": (agreement["fc_diff"] - agreement["fc_diff"].median()).abs().median(),
    }
