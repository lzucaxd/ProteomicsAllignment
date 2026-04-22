#!/usr/bin/env python3
"""
Legacy end-to-end benchmark (Python-chained steps).

Prefer the canonical shell orchestrator for slide-ready v2 output:
  ./scripts/benchmark/run_overnight_v2.sh
See docs/BENCHMARK_V2_AND_PRESENTATION.md and docs/END_TO_END_TECHNICAL_REPORT.md.
For MSstatsTMT (matrix build) vs limma (benchmark DA), see docs/INFERENCE_BASELINES.md.

This script runs, in order: preprocessing, native baselines, methods, benchmark, meeting exports.
Usage:
    python scripts/run_all.py
    python scripts/run_all.py --skip-preprocessing --skip-native
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

SCRIPTS_DIR = Path(__file__).resolve().parent


def run_step(name: str, script: str, extra_args: list[str] | None = None):
    """Run a pipeline step as a subprocess."""
    logger.info("\n" + "=" * 60)
    logger.info("  STEP: %s", name)
    logger.info("=" * 60)

    cmd = [sys.executable, str(SCRIPTS_DIR / script)]
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(cmd, cwd=str(SCRIPTS_DIR.parent))
    if result.returncode != 0:
        logger.error("Step '%s' failed with exit code %d", name, result.returncode)
        sys.exit(result.returncode)


def main():
    parser = argparse.ArgumentParser(description="End-to-end benchmark pipeline")
    parser.add_argument("--skip-preprocessing", action="store_true")
    parser.add_argument("--skip-native", action="store_true")
    parser.add_argument("--skip-methods", action="store_true")
    parser.add_argument("--skip-benchmark", action="store_true")
    parser.add_argument("--skip-exports", action="store_true")
    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("  CPTAC-CCLE HARMONIZATION BENCHMARK — FULL PIPELINE")
    logger.info("=" * 60)

    if not args.skip_preprocessing:
        run_step("Preprocessing", "run_preprocessing.py")

    if not args.skip_native:
        run_step("Native-domain baselines", "run_native_baselines.py")

    if not args.skip_methods:
        run_step("Method representations", "run_methods.py")

    if not args.skip_benchmark:
        run_step("Benchmark evaluation", "run_benchmark.py")

    if not args.skip_exports:
        run_step("Meeting exports", "run_meeting_exports.py")

    logger.info("\n" + "=" * 60)
    logger.info("  ALL STEPS COMPLETE")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
