#!/usr/bin/env python3
"""
Step 10: Build comparison_summary.csv from benchmark_results tree + cross_domain_metrics.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from harmonize.reporting.summary_tables import (  # noqa: E402
    annotate_metric_tiers,
    build_method_comparison_table,
)


def load_cross_domain_split(results: Path, method: str, task: str) -> dict:
    p = results / method / task / "representation_da" / "cross_domain_metrics.csv"
    out = {}
    if not p.is_file():
        return out
    df = pd.read_csv(p)
    for _, row in df.iterrows():
        lab = row.get("gene_set")
        if lab == "union":
            out["fc_correlation_union"] = row.get("fc_correlation")
            out["fc_same_dir_union"] = row.get("same_dir_fraction")
            out["n_genes_union"] = row.get("n_genes")
        elif lab == "intersection":
            out["fc_correlation_intersection"] = row.get("fc_correlation")
            out["fc_same_dir_intersection"] = row.get("same_dir_fraction")
            out["n_genes_intersection"] = row.get("n_genes")
    return out


def load_residual_neff(task_dir: Path) -> dict:
    out = {}
    for dom in ("cptac", "ccle"):
        f = task_dir / "calibration" / f"residual_dependence_{dom}.csv"
        if f.is_file():
            df = pd.read_csv(f)
            if len(df) > 0 and "effective_n" in df.columns:
                out[f"n_eff_{dom}"] = df["effective_n"].iloc[0]
            if len(df) > 0 and "mean_abs_residual_corr" in df.columns:
                out[f"residual_dependence_{dom}"] = df["mean_abs_residual_corr"].iloc[0]
    return out


def n_ccle_from_sample_meta(repo: Path, task: str) -> float:
    for rel in ("data/processed/union", "data/processed_union"):
        p = repo / rel / f"sample_meta_{task}.csv"
        if p.is_file():
            df = pd.read_csv(p)
            if "domain" in df.columns:
                return float((df["domain"].str.upper() == "CCLE").sum())
    return float("nan")


def load_disconnect_scores(results: Path) -> pd.DataFrame:
    p = results / "disconnect_scores.csv"
    if not p.is_file():
        return pd.DataFrame()
    return pd.read_csv(p)


def load_marker_counts(task_dir: Path, method: str, task: str) -> dict:
    cal = task_dir / "calibration"
    if not cal.is_dir():
        return {}
    out = {}
    for dom in ("cptac", "ccle"):
        for f in cal.glob(f"marker_sanity_summary_{method}_{dom}_{task}.csv"):
            df = pd.read_csv(f)
            if len(df) == 0:
                continue
            out[f"n_markers_tested_{dom}"] = int(df["n_markers_tested"].iloc[0])
            out[f"n_markers_correct_{dom}"] = int(df["n_correct"].iloc[0])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=REPO)
    ap.add_argument(
        "--results-root",
        type=Path,
        default=REPO / "reports" / "benchmark_master" / "benchmark_results",
    )
    args = ap.parse_args()
    results = args.results_root.resolve()

    base = build_method_comparison_table(results)
    if base.empty:
        print("No benchmark result dirs found", file=sys.stderr)
        sys.exit(1)

    extras = []
    for _, row in base.iterrows():
        method, task = row["method"], row["task"]
        task_dir = results / method / task
        d = load_cross_domain_split(results, method, task)
        d.update(load_residual_neff(task_dir))
        d.update(load_marker_counts(task_dir, method, task))
        # Prefer explicit union naming; keep legacy fc_correlation aligned with union
        if "fc_correlation_union" in d and pd.notna(d["fc_correlation_union"]):
            d["fc_correlation"] = d["fc_correlation_union"]
        if "fc_same_dir_union" in d and pd.notna(d["fc_same_dir_union"]):
            d["fc_same_dir_frac"] = d["fc_same_dir_union"]
        extras.append(d)

    extra_df = pd.DataFrame(extras)
    merged = pd.concat([base.reset_index(drop=True), extra_df], axis=1)

    merged["n_ccle_samples"] = merged["task"].map(lambda t: n_ccle_from_sample_meta(args.repo, str(t)))

    disc = load_disconnect_scores(results)
    if not disc.empty and all(c in disc.columns for c in ("method", "task", "disconnect_score")):
        merged = merged.merge(
            disc[["method", "task", "disconnect_score", "geom_improvement", "da_improvement"]],
            on=["method", "task"],
            how="left",
        )
    # Calibrated FC vs intersection ceiling
    if "fc_correlation_intersection" in merged.columns and "concordance_ceiling_fc_corr" in merged.columns:
        merged["calibrated_fc_corr_intersection"] = merged.apply(
            lambda r: r["fc_correlation_intersection"] / r["concordance_ceiling_fc_corr"]
            if pd.notna(r.get("fc_correlation_intersection"))
            and pd.notna(r.get("concordance_ceiling_fc_corr"))
            and r["concordance_ceiling_fc_corr"] > 0
            else np.nan,
            axis=1,
        )

    out_csv = results / "comparison_summary.csv"
    merged.to_csv(out_csv, index=False)
    print("Wrote", out_csv)

    tiered = annotate_metric_tiers(merged)
    tiered.to_csv(results / "comparison_summary_tiered.csv", index=False)


if __name__ == "__main__":
    main()
