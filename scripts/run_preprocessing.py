#!/usr/bin/env python3
"""
Preprocessing runner: loads raw data, builds shared spaces, saves benchmark-ready matrices.

Usage:
    python scripts/run_preprocessing.py [--config configs/preprocessing/default.yaml]
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[0] / ".." / "src"))

from harmonize.utils.config import load_config, load_task_config
from harmonize.utils.paths import ProjectPaths
from harmonize.utils.io import save_matrix, save_metadata
from harmonize.preprocessing.loaders import load_cptac_studies, load_ccle, load_ccle_sample_info
from harmonize.preprocessing.metadata import (
    average_cal120_ccle_columns,
    build_sample_meta,
    build_feature_meta,
)
from harmonize.preprocessing.shared_space import build_shared_space

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run preprocessing pipeline")
    parser.add_argument("--config", default="configs/preprocessing/default.yaml")
    parser.add_argument("--repo-root", default=".")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Override output_dir from config (e.g. data/processed_union)",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = ProjectPaths(args.repo_root)
    out_rel = args.output_dir or cfg.get("output_dir", "data/processed")
    out_dir = paths.ensure_dir(out_rel)

    logger.info("=" * 60)
    logger.info("  PREPROCESSING PIPELINE")
    logger.info("=" * 60)

    # ── Load data sources ───────────────────────────────────────────
    logger.info("\nStep 1: Loading data sources...")

    cptac_cfg = cfg.get("data_sources", {}).get("cptac", {})
    cptac_studies = {}
    for study_id, info in cptac_cfg.items():
        mat_path = info.get("gene_matrix")
        if mat_path and (paths.resolve(mat_path)).exists():
            cptac_studies[study_id] = {"tissue": info.get("tissue", "?"), "path": paths.resolve(mat_path)}

    cptac = load_cptac_studies(paths, cptac_studies)
    ccle = average_cal120_ccle_columns(load_ccle(paths))
    ccle_info = load_ccle_sample_info(paths)

    # ── Build shared spaces per task ────────────────────────────────
    shared_cfg = cfg.get("shared_space", {})
    filter_cfg = cfg.get("filtering", {})
    impute_cfg = cfg.get("imputation", {})

    tasks = ["breast_subtype", "breast_vs_lung"]

    for task_name in tasks:
        logger.info("\n" + "=" * 60)
        logger.info("  Task: %s", task_name)
        logger.info("=" * 60)

        task_cfg = load_task_config(task_name)

        # Select studies for this task
        if task_name == "breast_subtype":
            task_studies = {s: cptac[s] for s in task_cfg["cptac"]["studies"] if s in cptac}
        elif task_name == "breast_vs_lung":
            breast_s = task_cfg["cptac"].get("breast_study", "PDC000120")
            lung_s = task_cfg["cptac"].get("lung_study", "PDC000153")
            task_studies = {k: v for k, v in cptac.items() if k in (breast_s, lung_s)}
        else:
            task_studies = cptac

        # Build shared space
        cptac_shared, ccle_shared, stats = build_shared_space(
            task_studies,
            ccle,
            min_prevalence=filter_cfg.get("min_prevalence", 0.50),
            min_sd=filter_cfg.get("min_sd", 0.01),
            join_strategy=shared_cfg.get("join_strategy", "intersection"),
        )

        # Build metadata
        sample_meta = build_sample_meta(
            task_studies, ccle, task_cfg, ccle_info, repo_root=paths.root
        )
        sample_meta = sample_meta.drop_duplicates(subset=["sample_id"], keep="first")
        shared_genes = list(cptac_shared.index)
        feature_meta = build_feature_meta(shared_genes)

        # Combine into single matrix
        combined = cptac_shared.join(ccle_shared, how="inner")
        # Bidirectional filter: matrix ↔ metadata
        valid_sids = set(sample_meta["sample_id"])
        combined = combined[[c for c in combined.columns if c in valid_sids]]
        sample_meta = sample_meta[sample_meta["sample_id"].isin(combined.columns)]

        # Save
        prefix = task_name.replace(" ", "_")
        save_matrix(combined, out_dir / f"shared_gene_matrix_{prefix}.csv")
        save_metadata(sample_meta, out_dir / f"sample_meta_{prefix}.csv")
        save_metadata(feature_meta, out_dir / f"feature_meta_{prefix}.csv")

        logger.info("  Saved: shared_gene_matrix_%s.csv (%d genes x %d samples)",
                     prefix, len(combined), len(combined.columns))
        logger.info("  Saved: sample_meta_%s.csv (%d rows)", prefix, len(sample_meta))
        logger.info("  Stats: %s", stats)

    logger.info("\n" + "=" * 60)
    logger.info("  PREPROCESSING COMPLETE")
    logger.info("=" * 60)
    logger.info("Outputs in: %s", out_dir)


if __name__ == "__main__":
    main()
