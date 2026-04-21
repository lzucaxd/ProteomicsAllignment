"""Marker gene evaluation metrics."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def check_marker_availability(
    matrix: pd.DataFrame,
    markers: list[str],
) -> dict[str, bool]:
    """Check which markers are present in the matrix."""
    return {m: m in matrix.index for m in markers}


def extract_marker_summary(
    da_result: pd.DataFrame,
    markers: list[str],
    gene_col: str = "gene",
) -> pd.DataFrame:
    """Extract DA results for marker genes only."""
    if gene_col not in da_result.columns:
        for col in ("Protein", "GeneSymbol", "Gene"):
            if col in da_result.columns:
                gene_col = col
                break

    return da_result[da_result[gene_col].isin(markers)].copy()


def check_marker_directions(
    da_result: pd.DataFrame,
    expected_directions: dict[str, str],
    gene_col: str = "gene",
    fc_col: str = "logFC",
) -> pd.DataFrame:
    """
    Check whether marker fold-changes match expected directions.

    expected_directions: dict mapping gene -> "up_in_X" or "down_in_X"
    """
    records = []
    for gene, expected in expected_directions.items():
        rows = da_result[da_result[gene_col] == gene]
        if len(rows) == 0:
            records.append({"gene": gene, "expected": expected, "observed_fc": np.nan, "correct": None, "note": "not_in_da"})
            continue

        fc = rows.iloc[0][fc_col]
        if "up" in expected.lower():
            correct = fc > 0
        elif "down" in expected.lower():
            correct = fc < 0
        else:
            correct = None

        records.append({
            "gene": gene,
            "expected": expected,
            "observed_fc": fc,
            "correct": correct,
            "note": "",
        })

    return pd.DataFrame(records)


def compute_marker_concordance(
    da_cptac: pd.DataFrame,
    da_ccle: pd.DataFrame,
    markers: list[str],
    gene_col: str = "gene",
    fc_col: str = "logFC",
) -> dict:
    """Compute marker-level concordance between CPTAC and CCLE."""
    merged = pd.merge(
        da_cptac[[gene_col, fc_col]].rename(columns={fc_col: "fc_cptac"}),
        da_ccle[[gene_col, fc_col]].rename(columns={fc_col: "fc_ccle"}),
        on=gene_col,
    )
    marker_df = merged[merged[gene_col].isin(markers)]

    if len(marker_df) == 0:
        return {"n_markers_shared": 0}

    same_dir = (np.sign(marker_df["fc_cptac"]) == np.sign(marker_df["fc_ccle"])).sum()

    return {
        "n_markers_shared": len(marker_df),
        "n_same_direction": int(same_dir),
        "frac_same_direction": same_dir / len(marker_df),
        "marker_fc_correlation": marker_df[["fc_cptac", "fc_ccle"]].corr().iloc[0, 1],
    }
