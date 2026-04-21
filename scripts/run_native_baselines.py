#!/usr/bin/env python3
"""
Run native-domain baseline inference (Level 1).

Calls the existing R native_domain_da.R script for MSstatsTMT
and limma-based baselines on original, non-aligned data.

Usage:
    python scripts/run_native_baselines.py
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[0] / ".." / "src"))

from harmonize.utils.paths import ProjectPaths
from harmonize.utils.r_bridge import run_r_script

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run native-domain baselines")
    parser.add_argument("--repo-root", default=".")
    args = parser.parse_args()

    paths = ProjectPaths(args.repo_root)

    logger.info("=" * 60)
    logger.info("  NATIVE-DOMAIN BASELINES (Level 1)")
    logger.info("=" * 60)

    script = paths.r_benchmark_dir / "native_domain_da.R"
    if not script.exists():
        logger.error("Native DA script not found: %s", script)
        sys.exit(1)

    logger.info("Running native-domain DA...")
    run_r_script(script, cwd=paths.root)

    outdir = paths.resolve("reports/benchmark_master/native_domain_da")
    logger.info("Outputs in: %s", outdir)


if __name__ == "__main__":
    main()
