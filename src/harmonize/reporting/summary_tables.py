"""Consolidate benchmark results into comparison summary tables."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ── Metric tier definitions ─────────────────────────────────────────────────

METRIC_TIERS = {
    # High-confidence — directly interpretable after calibration
    "fc_correlation": "high",
    "fc_same_dir_frac": "high",
    "marker_sanity_cptac": "high",
    "marker_sanity_ccle": "high",
    "permutation_p_fc_corr": "high",
    # Informative — useful with context
    "struct_domain_r2_pc1": "informative",
    "struct_domain_r2_top5": "informative",
    "struct_condition_r2_pc1": "informative",
    "struct_classification_acc_domain": "informative",
    "struct_classification_acc_condition": "informative",
    "concordance_ceiling_fc_corr": "informative",
    "ccle_ceiling_fc_corr": "informative",
    "calibrated_fc_corr": "informative",
    "disconnect_score": "informative",
    "n_ccle_samples": "informative",
    "biology_destruction_retention": "informative",
    "biology_destruction_fc_shrinkage": "informative",
    # Supplementary — context-dependent
    "struct_silhouette_domain": "supplementary",
    "struct_silhouette_condition": "supplementary",
    "struct_knn_purity_domain": "supplementary",
    "struct_knn_purity_condition": "supplementary",
    "residual_mean_abs_corr": "supplementary",
    "residual_effective_n": "supplementary",
}


def annotate_metric_tiers(summary_df: pd.DataFrame) -> pd.DataFrame:
    """Add a 'tier' column based on METRIC_TIERS for each metric column.

    Returns a melted DataFrame with columns: method, task, metric, value, tier.
    """
    id_cols = ["method", "task"]
    metric_cols = [c for c in summary_df.columns if c not in id_cols]

    rows = []
    for _, row in summary_df.iterrows():
        for mc in metric_cols:
            val = row[mc]
            tier = METRIC_TIERS.get(mc, "unclassified")
            rows.append({
                "method": row.get("method"),
                "task": row.get("task"),
                "metric": mc,
                "value": val,
                "tier": tier,
            })

    return pd.DataFrame(rows)


def build_method_comparison_table(
    benchmark_results_dir: Path,
) -> pd.DataFrame:
    """
    Build a cross-method comparison table from individual benchmark results.

    Scans benchmark_results/<method>/<task>/ subdirectories for DA results,
    structure metrics, calibration outputs, and marker sanity files.
    """
    rows = []
    if not benchmark_results_dir.exists():
        return pd.DataFrame()

    for method_dir in sorted(benchmark_results_dir.iterdir()):
        if not method_dir.is_dir():
            continue
        for task_dir in sorted(method_dir.iterdir()):
            if not task_dir.is_dir():
                continue

            row = {"method": method_dir.name, "task": task_dir.name}

            # Structure metrics
            struct_file = task_dir / "structure" / "structure_summary.csv"
            if struct_file.exists():
                try:
                    df = pd.read_csv(struct_file)
                    if len(df) > 0:
                        for col in df.columns:
                            row[f"struct_{col}"] = df.iloc[0][col]
                except Exception:
                    pass

            # DA agreement
            agree_file = task_dir / "representation_da" / "fc_agreement.csv"
            if agree_file.exists():
                try:
                    df = pd.read_csv(agree_file)
                    row["fc_agreement_n"] = len(df)
                    if "same_direction" in df.columns:
                        row["fc_same_dir_frac"] = df["same_direction"].mean()
                    if "logFC_cptac" in df.columns and "logFC_ccle" in df.columns:
                        row["fc_correlation"] = df[["logFC_cptac", "logFC_ccle"]].corr().iloc[0, 1]
                except Exception:
                    pass

            # Calibration: permutation null
            perm_file = task_dir / "calibration" / "observed_vs_null_summary.csv"
            if perm_file.exists():
                try:
                    df = pd.read_csv(perm_file)
                    for _, prow in df.iterrows():
                        metric = prow.get("metric", "")
                        if metric == "fc_correlation":
                            row["permutation_p_fc_corr"] = prow.get("p_value")
                            row["permutation_z_fc_corr"] = prow.get("z_score")
                        elif metric == "same_direction_frac":
                            row["permutation_p_same_dir"] = prow.get("p_value")
                except Exception:
                    pass

            # Calibration: concordance ceiling (CPTAC primary; CCLE optional split-half)
            calib = task_dir / "calibration"
            ceil_cptac = calib / "ceiling_summary_cptac.csv"
            ceil_ccle = calib / "ceiling_summary_ccle.csv"
            ceil_legacy = calib / "ceiling_summary.csv"
            try:
                if ceil_cptac.exists():
                    df = pd.read_csv(ceil_cptac)
                elif ceil_legacy.exists():
                    df = pd.read_csv(ceil_legacy)
                else:
                    df = None
                if df is not None and "ceiling_fc_correlation" in df.columns:
                    row["concordance_ceiling_fc_corr"] = float(df["ceiling_fc_correlation"].iloc[0])
                    ceiling_val = row["concordance_ceiling_fc_corr"]
                    fc_corr = row.get("fc_correlation")
                    if fc_corr is not None and ceiling_val and ceiling_val > 0:
                        row["calibrated_fc_corr"] = fc_corr / ceiling_val
            except Exception:
                pass
            if ceil_ccle.exists():
                try:
                    df2 = pd.read_csv(ceil_ccle)
                    if "ceiling_fc_correlation" in df2.columns:
                        row["ccle_ceiling_fc_corr"] = float(df2["ceiling_fc_correlation"].iloc[0])
                except Exception:
                    pass

            # Calibration: marker sanity
            calib_dir = task_dir / "calibration"
            if calib_dir.exists():
                for domain in ["cptac", "ccle"]:
                    sanity_pattern = f"marker_sanity_summary_*_{domain}_*.csv"
                    for f in calib_dir.glob(sanity_pattern):
                        try:
                            df = pd.read_csv(f)
                            if "marker_sanity_rate" in df.columns:
                                row[f"marker_sanity_{domain}"] = df["marker_sanity_rate"].iloc[0]
                        except Exception:
                            pass

            # Calibration: biology destruction
            for f in (calib_dir.glob("destruction_summary_*.csv") if calib_dir.exists() else []):
                try:
                    df = pd.read_csv(f)
                    if len(df) > 0:
                        row["biology_destruction_retention"] = df.get(
                            "default_retention_rate", [np.nan]
                        ).iloc[0]
                        row["biology_destruction_fc_shrinkage"] = df.get(
                            "default_mean_fc_shrinkage", [np.nan]
                        ).iloc[0]
                except Exception:
                    pass

            # Calibration: residual dependence
            for f in (calib_dir.glob("residual_dependence_*.csv") if calib_dir.exists() else []):
                try:
                    df = pd.read_csv(f)
                    if len(df) > 0:
                        row["residual_mean_abs_corr"] = df.get(
                            "mean_abs_residual_corr", [np.nan]
                        ).iloc[0]
                        row["residual_effective_n"] = df.get(
                            "effective_n", [np.nan]
                        ).iloc[0]
                except Exception:
                    pass

            rows.append(row)

    return pd.DataFrame(rows)
