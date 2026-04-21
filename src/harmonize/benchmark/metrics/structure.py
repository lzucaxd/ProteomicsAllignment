"""Level 4: Structure and geometry evaluation metrics.

Quantifies domain mixing, biological structure preservation, and
residual domain effects in harmonized representations.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import silhouette_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import LabelEncoder

logger = logging.getLogger(__name__)


def compute_structure_metrics(
    matrix: pd.DataFrame,
    sample_meta: pd.DataFrame,
    config: dict[str, Any] | None = None,
    reference_matrix: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """
    Compute all structure / geometry metrics for a representation.

    Parameters
    ----------
    reference_matrix : optional raw/unharmonized matrix for fixed-basis PCA.
        If provided, PCA is fit on reference_matrix and the method matrix is
        projected onto that basis. This ensures all methods share the same
        PCA coordinate system for fair comparison.

    Metrics:
    - PCA variance decomposition (domain R², condition R² on top PCs)
    - kNN purity by domain and condition
    - Silhouette scores by domain and condition
    - Domain / condition classification accuracy (logistic regression on PCs)
    """
    config = config or {}
    n_pcs = config.get("n_pcs", 20)
    knn_k = config.get("knn_k", 15)

    # Align matrix columns with metadata
    meta = sample_meta.copy()
    valid_cols = [c for c in matrix.columns if c in set(meta["sample_id"])]
    if len(valid_cols) < 10:
        return {"error": f"Too few matching samples ({len(valid_cols)})"}

    meta = meta[meta["sample_id"].isin(valid_cols)].set_index("sample_id")
    X = matrix[valid_cols].T
    X = X.fillna(0)
    meta = meta.loc[valid_cols]

    domain = meta["domain"].values if "domain" in meta.columns else None
    condition = meta["condition"].values if "condition" in meta.columns else None

    results = {"n_samples": len(valid_cols), "n_genes": len(matrix)}

    # PCA — fixed-basis or per-method
    actual_pcs = min(n_pcs, min(X.shape) - 1)
    if actual_pcs < 2:
        return {**results, "error": "Not enough dimensions for PCA"}

    pca = PCA(n_components=actual_pcs)

    if reference_matrix is not None:
        # Fixed-basis PCA: fit on reference (raw) data, project method data
        shared_genes = sorted(set(matrix.index) & set(reference_matrix.index))
        ref_cols = [c for c in reference_matrix.columns if c in set(meta.index)]
        if len(shared_genes) < 50 or len(ref_cols) < 10:
            logger.warning("Fixed-basis PCA: insufficient overlap, falling back to per-method")
            Z = pca.fit_transform(X.values)
            results["pca_basis"] = "per_method"
        else:
            X_ref = reference_matrix.loc[shared_genes, ref_cols].T.fillna(0)
            X_method = matrix.loc[shared_genes, valid_cols].T.fillna(0)
            pca.fit(X_ref.values)
            Z = pca.transform(X_method.values)
            results["pca_basis"] = "fixed"
    else:
        Z = pca.fit_transform(X.values)
        results["pca_basis"] = "per_method"

    results["pca_variance_explained"] = pca.explained_variance_ratio_.tolist()

    # R² on PCs (ANOVA-style)
    if domain is not None:
        results["domain_r2_pc1"] = _r2_from_labels(Z[:, 0], domain)
        results["domain_r2_pc2"] = _r2_from_labels(Z[:, 1], domain)
        results["domain_r2_top5"] = _multivariate_r2(Z[:, :5], domain)

    if condition is not None:
        results["condition_r2_pc1"] = _r2_from_labels(Z[:, 0], condition)
        results["condition_r2_pc2"] = _r2_from_labels(Z[:, 1], condition)
        results["condition_r2_top5"] = _multivariate_r2(Z[:, :5], condition)

    # Silhouette scores
    if domain is not None and len(set(domain)) >= 2:
        try:
            results["silhouette_domain"] = float(silhouette_score(Z[:, :actual_pcs], domain))
        except Exception:
            results["silhouette_domain"] = float("nan")

    if condition is not None and len(set(condition)) >= 2:
        try:
            results["silhouette_condition"] = float(silhouette_score(Z[:, :actual_pcs], condition))
        except Exception:
            results["silhouette_condition"] = float("nan")

    # kNN purity
    if domain is not None and len(set(domain)) >= 2:
        results["knn_purity_domain"] = _knn_purity(Z[:, :actual_pcs], domain, knn_k)

    if condition is not None and len(set(condition)) >= 2:
        results["knn_purity_condition"] = _knn_purity(Z[:, :actual_pcs], condition, knn_k)

    # Classification accuracy (logistic regression on PCs)
    if domain is not None and len(set(domain)) >= 2:
        results["classification_acc_domain"] = _logreg_accuracy(Z[:, :actual_pcs], domain)

    if condition is not None and len(set(condition)) >= 2:
        results["classification_acc_condition"] = _logreg_accuracy(Z[:, :actual_pcs], condition)

    return results


def _r2_from_labels(values: np.ndarray, labels: np.ndarray) -> float:
    """One-way ANOVA R² (SS_between / SS_total)."""
    overall_mean = np.mean(values)
    ss_total = np.sum((values - overall_mean) ** 2)
    if ss_total == 0:
        return 0.0

    unique = np.unique(labels)
    ss_between = 0.0
    for lab in unique:
        mask = labels == lab
        group_mean = np.mean(values[mask])
        ss_between += mask.sum() * (group_mean - overall_mean) ** 2

    return float(ss_between / ss_total)


def _multivariate_r2(Z: np.ndarray, labels: np.ndarray) -> float:
    """Average R² across multiple PCs."""
    r2s = [_r2_from_labels(Z[:, i], labels) for i in range(Z.shape[1])]
    return float(np.mean(r2s))


def _knn_purity(Z: np.ndarray, labels: np.ndarray, k: int = 15) -> float:
    """Fraction of k-nearest neighbors sharing the same label (average over all points)."""
    from sklearn.neighbors import NearestNeighbors

    k = min(k, len(Z) - 1)
    if k < 1:
        return float("nan")

    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")
    nn.fit(Z)
    distances, indices = nn.kneighbors(Z)

    purities = []
    for i in range(len(Z)):
        neighbors = indices[i, 1:]  # exclude self
        same = np.sum(labels[neighbors] == labels[i])
        purities.append(same / k)

    return float(np.mean(purities))


def _logreg_accuracy(Z: np.ndarray, labels: np.ndarray) -> float:
    """Leave-one-out-ish logistic regression classification accuracy."""
    le = LabelEncoder()
    y = le.fit_transform(labels)

    if len(np.unique(y)) < 2:
        return float("nan")

    try:
        clf = LogisticRegression(max_iter=1000, solver="lbfgs", multi_class="auto")
        clf.fit(Z, y)
        return float(clf.score(Z, y))
    except Exception:
        return float("nan")
