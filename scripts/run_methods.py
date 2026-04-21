#!/usr/bin/env python3
"""
Method runner: generate harmonized representations for the benchmark.

Usage:
    python scripts/run_methods.py [--methods raw bridge_aware celligner]
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[0] / ".." / "src"))

from harmonize.utils.config import load_config, load_method_config
from harmonize.utils.paths import ProjectPaths
from harmonize.methods.registry import get_method, list_methods

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run harmonization methods")
    parser.add_argument("--methods", nargs="+", default=None,
                        help="Methods to run (default: all from benchmark config)")
    parser.add_argument("--benchmark-config", default="configs/benchmark/default.yaml")
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip methods whose output already exists")
    args = parser.parse_args()

    paths = ProjectPaths(args.repo_root)
    bench_cfg = load_config(args.benchmark_config)

    methods_to_run = args.methods or bench_cfg.get("methods", list_methods())

    logger.info("=" * 60)
    logger.info("  METHOD RUNNER")
    logger.info("=" * 60)
    logger.info("Methods: %s", ", ".join(methods_to_run))

    for method_name in methods_to_run:
        logger.info("\n" + "-" * 50)
        logger.info("  Method: %s", method_name)
        logger.info("-" * 50)

        try:
            method_cfg = load_method_config(method_name)
        except FileNotFoundError:
            logger.warning("No config found for method '%s', skipping", method_name)
            continue

        outdir = paths.resolve(method_cfg.get("output_dir", f"reports/benchmark_master/methods/{method_name}"))

        if args.skip_existing and (outdir / "transformed_matrix.csv").exists():
            logger.info("  Output exists, skipping (--skip-existing)")
            continue

        method = get_method(method_name, paths)
        logger.info("  Output dir: %s", outdir)
        logger.info("  Instantiated: %s (%s)", method.display_name, type(method).__name__)

        # Methods that need re-running call method.run(...)
        # For now, just verify the method can be loaded
        logger.info("  Method '%s' registered and ready.", method_name)

    logger.info("\n" + "=" * 60)
    logger.info("  METHOD RUNNER COMPLETE")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
