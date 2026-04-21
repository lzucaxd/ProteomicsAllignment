"""Build shared feature space across CPTAC and CCLE."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from harmonize.utils.io import intersect_genes, union_genes, filter_by_prevalence, remove_near_constant

logger = logging.getLogger(__name__)


def build_shared_space(
    cptac_matrices: dict[str, pd.DataFrame],
    ccle_matrix: pd.DataFrame,
    min_prevalence: float = 0.50,
    min_sd: float = 0.01,
    join_strategy: str = "intersection",
    impute: bool = False,
    impute_strategy: str = "within_domain_gene_median",
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Build a shared gene-level matrix from CPTAC + CCLE.

    Parameters
    ----------
    join_strategy : "intersection" or "union"
        "intersection" keeps only genes present in both domains (conservative,
        no missing data introduced at this step).
        "union" keeps all genes from either domain (maximizes feature coverage
        but introduces NaN for genes missing in one domain — pair with
        imputation or downstream filtering).

    Steps:
    1. Union CPTAC studies (outer join across studies)
    2. Intersect or union with CCLE genes
    3. Prevalence filtering within each domain
    4. Remove near-constant genes
    5. Optional imputation

    Returns (cptac_shared, ccle_shared, stats_dict).
    """
    if join_strategy not in ("intersection", "union"):
        raise ValueError(f"join_strategy must be 'intersection' or 'union', got '{join_strategy}'")

    # Union CPTAC studies
    cptac_dfs = list(cptac_matrices.values())
    if len(cptac_dfs) == 1:
        cptac_union = cptac_dfs[0]
    else:
        cptac_union = pd.concat(cptac_dfs, axis=1, join="outer")

    # Join gene sets across domains
    if join_strategy == "intersection":
        shared = intersect_genes(cptac_union, ccle_matrix)
        logger.info("Gene join (intersection): %d genes", len(shared))
    else:
        shared = union_genes(cptac_union, ccle_matrix)
        n_cptac_only = len(set(cptac_union.index) - set(ccle_matrix.index))
        n_ccle_only = len(set(ccle_matrix.index) - set(cptac_union.index))
        logger.info("Gene join (union): %d genes (%d CPTAC-only, %d CCLE-only)",
                     len(shared), n_cptac_only, n_ccle_only)

    cptac_shared = cptac_union.reindex(shared)
    ccle_shared = ccle_matrix.reindex(shared)

    # Prevalence filtering per domain
    cptac_prev, _cptac_dropped = filter_by_prevalence(cptac_shared, min_prevalence)
    ccle_prev, _ccle_dropped = filter_by_prevalence(ccle_shared, min_prevalence)

    if join_strategy == "intersection":
        genes_after_prev = sorted(set(cptac_prev.index) & set(ccle_prev.index))
    else:
        # Union: keep genes that pass prevalence in at least one domain (union-specific
        # genes remain as NaN in the other domain).
        genes_after_prev = sorted(set(cptac_prev.index) | set(ccle_prev.index))
    logger.info(
        "After prevalence filter (%.0f%%, %s): %d genes",
        min_prevalence * 100,
        join_strategy,
        len(genes_after_prev),
    )

    cptac_shared = cptac_shared.reindex(genes_after_prev)
    ccle_shared = ccle_shared.reindex(genes_after_prev)

    # Remove near-constant
    combined = pd.concat([cptac_shared, ccle_shared], axis=1)
    combined = remove_near_constant(combined, min_sd)
    genes_final = list(combined.index)
    logger.info("After SD filter (>%.3f): %d genes", min_sd, len(genes_final))

    cptac_shared = cptac_shared.loc[genes_final]
    ccle_shared = ccle_shared.loc[genes_final]

    # Optional imputation
    if impute and impute_strategy == "within_domain_gene_median":
        cptac_shared = _impute_median(cptac_shared, "CPTAC")
        ccle_shared = _impute_median(ccle_shared, "CCLE")

    stats = {
        "join_strategy": join_strategy,
        "n_joined_raw": len(shared),
        "n_after_prevalence": len(genes_after_prev),
        "n_final": len(genes_final),
        "n_cptac_samples": len(cptac_shared.columns),
        "n_ccle_samples": len(ccle_shared.columns),
    }

    return cptac_shared, ccle_shared, stats


def _impute_median(df: pd.DataFrame, label: str) -> pd.DataFrame:
    """Within-domain gene-median imputation."""
    n_missing = df.isna().sum().sum()
    if n_missing > 0:
        medians = df.median(axis=1)
        df = df.T.fillna(medians).T
        logger.info("Imputed %d missing values in %s (gene median)", n_missing, label)
    return df
