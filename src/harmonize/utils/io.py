"""I/O utilities for gene matrices and metadata tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_gene_matrix(path: str | Path) -> pd.DataFrame:
    """
    Load a gene matrix CSV (genes x samples).

    Returns a DataFrame with gene symbols as the index and sample IDs as columns.
    Drops non-numeric ID columns (GeneSymbol, UniProtID, Gene).
    """
    df = pd.read_csv(path, index_col=0)
    id_cols = [c for c in df.columns if c in ("UniProtID", "Gene")]
    if id_cols:
        df = df.drop(columns=id_cols)
    df.index = df.index.astype(str)
    df.index.name = "GeneSymbol"
    return df.apply(pd.to_numeric, errors="coerce")


def load_metadata(path: str | Path) -> pd.DataFrame:
    """Load a sample or feature metadata CSV."""
    return pd.read_csv(path)


def save_matrix(df: pd.DataFrame, path: str | Path) -> None:
    """Save a gene matrix to CSV with index."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path)


def save_metadata(df: pd.DataFrame, path: str | Path) -> None:
    """Save metadata DataFrame to CSV without index."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def intersect_genes(*matrices: pd.DataFrame) -> list[str]:
    """Return sorted list of gene symbols present in all matrices."""
    if not matrices:
        return []
    common = set(matrices[0].index)
    for m in matrices[1:]:
        common &= set(m.index)
    return sorted(common)


def union_genes(*matrices: pd.DataFrame) -> list[str]:
    """Return sorted list of gene symbols present in any matrix."""
    if not matrices:
        return []
    all_genes: set[str] = set()
    for m in matrices:
        all_genes |= set(m.index)
    return sorted(all_genes)


def filter_by_prevalence(
    df: pd.DataFrame, min_frac: float = 0.50
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Keep genes observed in at least `min_frac` of samples.

    Returns (filtered_df, dropped_df).
    """
    obs_frac = df.notna().mean(axis=1)
    mask = obs_frac >= min_frac
    return df.loc[mask], df.loc[~mask]


def remove_near_constant(df: pd.DataFrame, min_sd: float = 0.01) -> pd.DataFrame:
    """Remove genes with standard deviation below threshold."""
    sds = df.std(axis=1, skipna=True)
    return df.loc[sds >= min_sd]
