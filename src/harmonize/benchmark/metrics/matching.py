"""Level 5: Tumor-cell line matching / retrieval evaluation.

Evaluates whether harmonization improves the ability to match tumors
to their most biologically similar cell lines.

STATUS: Scaffold with partial implementation.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_matching_metrics(
    matrix: pd.DataFrame,
    sample_meta: pd.DataFrame,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Compute tumor-cell line matching quality metrics.

    For each CPTAC tumor sample, find its k nearest CCLE neighbors
    and evaluate whether they share the expected lineage/subtype.
    """
    config = config or {}
    top_k_values = config.get("top_k", [1, 3, 5, 10])
    distance_metric = config.get("distance_metric", "correlation")

    meta = sample_meta.copy()
    valid_cols = [c for c in matrix.columns if c in set(meta["sample_id"])]
    meta = meta[meta["sample_id"].isin(valid_cols)].set_index("sample_id")

    cptac_ids = meta[meta["domain"] == "CPTAC"].index.tolist()
    ccle_ids = meta[meta["domain"] == "CCLE"].index.tolist()

    if len(cptac_ids) < 5 or len(ccle_ids) < 5:
        return {"error": "Too few samples for matching", "n_cptac": len(cptac_ids), "n_ccle": len(ccle_ids)}

    cptac_ids = [c for c in cptac_ids if c in matrix.columns]
    ccle_ids = [c for c in ccle_ids if c in matrix.columns]

    X_cptac = matrix[cptac_ids].T.fillna(0).values
    X_ccle = matrix[ccle_ids].T.fillna(0).values

    # Compute pairwise distances
    dist_matrix = _pairwise_distances(X_cptac, X_ccle, metric=distance_metric)

    # For each tumor, rank cell lines by distance
    results = {"n_cptac": len(cptac_ids), "n_ccle": len(ccle_ids)}

    if "condition" in meta.columns:
        cptac_conditions = meta.loc[cptac_ids, "condition"].values
        ccle_conditions = meta.loc[ccle_ids, "condition"].values

        for k in top_k_values:
            if k > len(ccle_ids):
                continue
            match_rate = _top_k_match_rate(dist_matrix, cptac_conditions, ccle_conditions, k)
            results[f"top_{k}_same_condition_rate"] = match_rate

    return results


def _pairwise_distances(
    X: np.ndarray,
    Y: np.ndarray,
    metric: str = "correlation",
) -> np.ndarray:
    """Compute pairwise distance matrix between rows of X and Y."""
    from scipy.spatial.distance import cdist

    if metric == "correlation":
        X_centered = X - X.mean(axis=1, keepdims=True)
        Y_centered = Y - Y.mean(axis=1, keepdims=True)
        X_std = np.linalg.norm(X_centered, axis=1, keepdims=True)
        Y_std = np.linalg.norm(Y_centered, axis=1, keepdims=True)
        X_std[X_std == 0] = 1
        Y_std[Y_std == 0] = 1
        X_norm = X_centered / X_std
        Y_norm = Y_centered / Y_std
        similarity = X_norm @ Y_norm.T
        return 1 - similarity

    return cdist(X, Y, metric=metric)


def _top_k_match_rate(
    dist_matrix: np.ndarray,
    query_labels: np.ndarray,
    reference_labels: np.ndarray,
    k: int,
) -> float:
    """Fraction of queries whose top-k neighbors include at least one same-label match."""
    n_queries = dist_matrix.shape[0]
    matches = 0

    for i in range(n_queries):
        ranked_indices = np.argsort(dist_matrix[i])[:k]
        neighbor_labels = reference_labels[ranked_indices]
        if query_labels[i] in neighbor_labels:
            matches += 1

    return matches / n_queries


def compute_nn_consistency(
    matrices: dict[str, pd.DataFrame],
    sample_meta: pd.DataFrame,
    k: int = 5,
) -> pd.DataFrame:
    """
    Compare nearest-neighbor assignments across methods.

    For each tumor sample, check whether its top-k CCLE neighbors
    are consistent across different method representations.

    Returns DataFrame with overlap statistics.
    """
    meta = sample_meta.copy()
    if "domain" not in meta.columns:
        return pd.DataFrame()

    cptac_ids = meta[meta["domain"] == "CPTAC"]["sample_id"].tolist()
    ccle_ids = meta[meta["domain"] == "CCLE"]["sample_id"].tolist()

    if len(cptac_ids) < 2 or len(ccle_ids) < 2:
        return pd.DataFrame()

    nn_per_method = {}
    for method_name, matrix in matrices.items():
        valid_cptac = [c for c in cptac_ids if c in matrix.columns]
        valid_ccle = [c for c in ccle_ids if c in matrix.columns]
        if len(valid_cptac) < 2 or len(valid_ccle) < 2:
            continue

        X_c = matrix[valid_cptac].T.fillna(0).values
        X_e = matrix[valid_ccle].T.fillna(0).values
        dists = _pairwise_distances(X_c, X_e)

        nn_sets = {}
        for i, sid in enumerate(valid_cptac):
            top_k_idx = np.argsort(dists[i])[:k]
            nn_sets[sid] = set(np.array(valid_ccle)[top_k_idx])
        nn_per_method[method_name] = nn_sets

    # Compute pairwise overlap between methods
    method_names = list(nn_per_method.keys())
    records = []
    for i in range(len(method_names)):
        for j in range(i + 1, len(method_names)):
            m1, m2 = method_names[i], method_names[j]
            shared_sids = set(nn_per_method[m1].keys()) & set(nn_per_method[m2].keys())
            overlaps = []
            for sid in shared_sids:
                s1 = nn_per_method[m1][sid]
                s2 = nn_per_method[m2][sid]
                overlap = len(s1 & s2) / k
                overlaps.append(overlap)
            records.append({
                "method_1": m1,
                "method_2": m2,
                "n_tumors": len(shared_sids),
                "mean_nn_overlap": np.mean(overlaps) if overlaps else float("nan"),
            })

    return pd.DataFrame(records)
