"""Imputation strategies for benchmark matrices."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def impute_within_domain_median(df: pd.DataFrame, label: str = "") -> pd.DataFrame:
    """Replace NaN with per-gene median (computed within this matrix only)."""
    n_before = df.isna().sum().sum()
    if n_before == 0:
        return df
    medians = df.median(axis=1)
    result = df.T.fillna(medians).T
    logger.info("Imputed %d values in %s (within-domain gene median)", n_before, label)
    return result


def zscore_standardize(df: pd.DataFrame) -> pd.DataFrame:
    """Row-wise z-score standardization (per gene)."""
    means = df.mean(axis=1)
    sds = df.std(axis=1)
    sds = sds.replace(0, 1.0)
    return df.sub(means, axis=0).div(sds, axis=0)
