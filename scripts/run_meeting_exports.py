#!/usr/bin/env python3
"""
Export meeting-ready materials from benchmark outputs.

Usage:
    python scripts/run_meeting_exports.py [--output reports/benchmark_master/meeting/export]
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[0] / ".." / "src"))

from harmonize.utils.paths import ProjectPaths
from harmonize.reporting.meeting_exports import export_meeting_materials
from harmonize.reporting.summary_tables import build_method_comparison_table

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Export meeting materials")
    parser.add_argument("--output", default="reports/benchmark_master/meeting/export")
    parser.add_argument("--repo-root", default=".")
    args = parser.parse_args()

    paths = ProjectPaths(args.repo_root)
    outdir = paths.resolve(args.output)

    logger.info("=" * 60)
    logger.info("  MEETING EXPORT")
    logger.info("=" * 60)

    exported = export_meeting_materials(paths.benchmark_root, outdir)

    # Build comparison table
    results_dir = paths.resolve("reports/benchmark_master/benchmark_results")
    comparison = build_method_comparison_table(results_dir)
    if len(comparison) > 0:
        comp_path = outdir / "tables" / "method_comparison.csv"
        comp_path.parent.mkdir(parents=True, exist_ok=True)
        comparison.to_csv(comp_path, index=False)
        logger.info("Comparison table: %s (%d rows)", comp_path, len(comparison))

    logger.info("\nExported %d files to %s", len(exported), outdir)


if __name__ == "__main__":
    main()
