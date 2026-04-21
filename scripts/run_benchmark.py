#!/usr/bin/env python3
"""
Benchmark runner: evaluate all methods on all tasks with calibration.

Usage:
    python scripts/run_benchmark.py
    python scripts/run_benchmark.py --methods raw bridge_shift --tasks breast_subtype
    python scripts/run_benchmark.py --skip-calibration

Union gene space (after: python scripts/run_preprocessing.py --config configs/preprocessing/union.yaml):
    python scripts/run_benchmark.py --processed-dir data/processed_union \\
        --benchmark-results-dir reports/benchmark_master/benchmark_results_union
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[0] / ".." / "src"))

from harmonize.utils.config import load_config, load_task_config, load_method_config
from harmonize.utils.paths import ProjectPaths
from harmonize.utils.io import load_gene_matrix
from harmonize.utils.r_bridge import run_r_script
from harmonize.benchmark.tasks import TaskDefinition
from harmonize.benchmark.runner import run_benchmark
from harmonize.methods.base import MethodResult
from harmonize.methods.bridge_aware import BridgeAwareMethod
from harmonize.methods.celligner import CellignerMethod

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


# =============================================================================
# Method loading
# =============================================================================

def _load_processed_matrix_if_present(
    paths: ProjectPaths, processed_dir: Path, task_name: str, sample_meta,
) -> "tuple[pd.DataFrame | None, list[str] | None]":
    """If preprocessing outputs exist, return (genes x samples matrix, gene list)."""
    import pandas as pd
    from harmonize.utils.io import load_gene_matrix

    mat_path = processed_dir / f"shared_gene_matrix_{task_name}.csv"
    if not mat_path.exists():
        return None, None
    combined = load_gene_matrix(mat_path)
    valid = set(sample_meta["sample_id"])
    cols = [c for c in combined.columns if c in valid]
    combined = combined[cols]
    return combined, combined.index.tolist()


def _align_matrix_genes(result: MethodResult, gene_index: list[str]) -> MethodResult:
    """Reindex method matrix rows to match raw/preprocessed gene order (NaN for missing)."""
    import pandas as pd

    result.matrix = result.matrix.reindex(gene_index)
    result.feature_meta = pd.DataFrame({"gene": gene_index, "included": True})
    return result


def _load_method_result(
    method_key: str,
    task: TaskDefinition,
    paths: ProjectPaths,
    sample_meta,
    processed_dir: Path,
    align_gene_index: list[str] | None = None,
) -> MethodResult | None:
    """Load a pre-computed method result from existing outputs."""
    import pandas as pd

    try:
        if method_key == "raw":
            from harmonize.preprocessing.loaders import load_cptac_studies, load_ccle
            from harmonize.utils.io import intersect_genes

            combined, genes = _load_processed_matrix_if_present(
                paths, processed_dir, task.name, sample_meta,
            )
            if combined is not None:
                logger.info(
                    "  Raw: loaded preprocessed matrix (%d genes) from %s",
                    len(genes),
                    processed_dir,
                )
                return MethodResult(
                    matrix=combined,
                    sample_meta=sample_meta,
                    feature_meta=pd.DataFrame({"gene": genes, "included": True}),
                    method_name="raw",
                    display_name="Raw",
                )

            cptac = load_cptac_studies(paths)
            ccle = load_ccle(paths)
            if task.name == "breast_subtype":
                study = task.raw_config["cptac"]["studies"][0]
                mats = [cptac[study]] if study in cptac else list(cptac.values())
            else:
                mats = list(cptac.values())
            shared = intersect_genes(*mats, ccle)
            combined = pd.concat([m.loc[shared] for m in mats] + [ccle.loc[shared]], axis=1)
            valid = set(sample_meta["sample_id"])
            combined = combined[[c for c in combined.columns if c in valid]]
            return MethodResult(
                matrix=combined,
                sample_meta=sample_meta,
                feature_meta=pd.DataFrame({"gene": shared, "included": True}),
                method_name="raw",
                display_name="Raw",
            )

        elif method_key in ("bridge_shift", "bridge_scale"):
            bridge = BridgeAwareMethod(paths)
            mode = "shift_only" if "shift" in method_key else "shift_and_scale"
            result = bridge.load_existing(mode)
            if result:
                result.sample_meta = sample_meta
                if align_gene_index is not None:
                    result = _align_matrix_genes(result, align_gene_index)
            return result

        elif method_key == "celligner":
            cell = CellignerMethod(paths)
            result = cell.load_existing(sample_meta)
            if align_gene_index is not None:
                result = _align_matrix_genes(result, align_gene_index)
            return result

    except Exception as e:
        logger.error("Failed to load %s: %s", method_key, e)
        return None


def _load_raw_reference(
    task: TaskDefinition,
    paths: ProjectPaths,
    sample_meta,
    processed_dir: Path,
) -> "pd.DataFrame | None":
    """Load raw shared matrix for fixed-basis PCA reference."""
    try:
        result = _load_method_result("raw", task, paths, sample_meta, processed_dir)
        if result is not None:
            return result.matrix
    except Exception as e:
        logger.warning("Could not load raw reference matrix: %s", e)
    return None


# =============================================================================
# Calibration helpers
# =============================================================================

def _parse_contrast(task: TaskDefinition) -> tuple[str, str]:
    """Parse contrast string into (contrast_a, contrast_b)."""
    parts = task.contrast.split("_vs_")
    if len(parts) == 2:
        return parts[1], parts[0]  # contrast_a, contrast_b
    return "GroupA", "GroupB"


def _build_expected_directions(task: TaskDefinition) -> tuple[str, str] | None:
    """Build marker/expected-sign CSV strings from task config."""
    if not task.expected_directions:
        return None

    contrast_a, contrast_b = _parse_contrast(task)
    markers = []
    signs = []

    for gene, direction in task.expected_directions.items():
        if direction == "variable":
            continue
        markers.append(gene)
        # Positive logFC = higher in contrast_b (the first part of X_vs_Y)
        if "up_in" in direction:
            up_label = direction.replace("up_in_", "").lower()
            if up_label == contrast_b.lower():
                signs.append("1")
            elif up_label == contrast_a.lower():
                signs.append("-1")
            else:
                signs.append("1")
        else:
            signs.append("1")

    if not markers:
        return None
    return ",".join(markers), ",".join(signs)


def _run_contrast_validation(
    task: TaskDefinition, paths: ProjectPaths, processed_dir: Path,
) -> bool:
    """Phase A: Run contrast validation on raw subtype data (task breast_subtype only)."""
    if task.name != "breast_subtype":
        return True

    logger.info("  Phase A: Contrast direction validation")

    diag_script = paths.r_benchmark_dir / "diagnose_subtype_sign.R"
    if not diag_script.exists():
        logger.warning("  diagnose_subtype_sign.R not found, skipping validation")
        return True

    processed_meta = processed_dir / "sample_meta_breast_subtype.csv"
    processed_mat = processed_dir / "shared_gene_matrix_breast_subtype.csv"

    if not processed_meta.exists() or not processed_mat.exists():
        logger.warning("  Preprocessed subtype data not found under %s, skipping validation", processed_dir)
        return True

    try:
        run_r_script(
            diag_script,
            args=[
                "--repo-root", str(paths.root),
                "--matrix", str(processed_mat),
                "--meta", str(processed_meta),
            ],
            cwd=paths.root,
        )
    except Exception as e:
        logger.error("  Contrast validation failed: %s", e)
        return False

    # Check the output
    import pandas as pd
    diag_file = paths.resolve("reports/benchmark_master/diagnostics/subtype_sign_diagnostic_summary.csv")
    if diag_file.exists():
        diag = pd.read_csv(diag_file)
        if "likely_flipped" in diag.columns and diag["likely_flipped"].any():
            logger.error("  CONTRAST LIKELY FLIPPED — check diagnostics before proceeding")
            return False
        logger.info("  Contrast validation PASSED for all domains")
    return True


def _run_calibration_for_method(
    method_key: str,
    task: TaskDefinition,
    paths: ProjectPaths,
    outdir: Path,
    sample_meta,
    processed_dir: Path,
    raw_outdir: Path | None = None,
):
    """Phase B: Run per-method calibration modules (permutation, ceiling, etc.)."""
    import pandas as pd

    calib_dir = outdir / "calibration"
    calib_dir.mkdir(parents=True, exist_ok=True)

    contrast_a, contrast_b = _parse_contrast(task)
    marker_info = _build_expected_directions(task)

    # Locate the method's DA results for input
    da_dir = outdir / "representation_da"

    # Write shared matrix + meta for R calibration scripts
    mat_path = processed_dir / f"shared_gene_matrix_{task.name}.csv"
    meta_path = processed_dir / f"sample_meta_{task.name}.csv"
    if not mat_path.exists() or not meta_path.exists():
        logger.warning("  Calibration skipped: preprocessed data not found for %s", task.name)
        return

    calibration_r_dir = paths.root / "src" / "harmonize" / "benchmark" / "calibration"

    # ── Permutation null ─────────────────────────────────────────────────
    perm_script = calibration_r_dir / "permutation_null.R"
    if perm_script.exists():
        logger.info("    Permutation null calibration (%s x %s)", method_key, task.name)
        perm_args = [
            "--matrix", str(da_dir / "_input_matrix.csv") if (da_dir / "_input_matrix.csv").exists() else str(mat_path),
            "--meta", str(da_dir / "_input_meta.csv") if (da_dir / "_input_meta.csv").exists() else str(meta_path),
            "--contrast-a", contrast_a,
            "--contrast-b", contrast_b,
            "--n-perm", "200",
            "--seed", "42",
            "--outdir", str(calib_dir),
        ]
        if marker_info:
            perm_args.extend(["--markers", marker_info[0], "--expected-signs", marker_info[1]])
        try:
            run_r_script(perm_script, args=perm_args, cwd=paths.root, timeout=7200)
        except Exception as e:
            logger.error("    Permutation null failed: %s", e)

    # ── Concordance ceiling (per domain) ────────────────────────────────
    ceil_script = calibration_r_dir / "concordance_ceiling.R"
    if ceil_script.exists():
        for domain in ["CPTAC", "CCLE"]:
            logger.info("    Concordance ceiling: %s", domain)
            try:
                run_r_script(ceil_script, args=[
                    "--matrix", str(mat_path),
                    "--meta", str(meta_path),
                    "--domain", domain,
                    "--contrast-a", contrast_a,
                    "--contrast-b", contrast_b,
                    "--n-splits", "100",
                    "--seed", "42",
                    "--outdir", str(calib_dir),
                ], cwd=paths.root, timeout=3600)
            except Exception as e:
                logger.error("    Concordance ceiling failed for %s: %s", domain, e)

    # ── Residual dependence (per domain) ────────────────────────────────
    resid_script = calibration_r_dir / "residual_dependence.R"
    if resid_script.exists():
        for domain in ["CPTAC", "CCLE"]:
            logger.info("    Residual dependence: %s", domain)
            try:
                run_r_script(resid_script, args=[
                    "--matrix", str(mat_path),
                    "--meta", str(meta_path),
                    "--contrast-a", contrast_a,
                    "--contrast-b", contrast_b,
                    "--domain", domain,
                    "--outdir", str(calib_dir),
                ], cwd=paths.root)
            except Exception as e:
                logger.error("    Residual dependence failed for %s: %s", domain, e)

    # ── Biology destruction ─────────────────────────────────────────────
    destr_script = calibration_r_dir / "biology_destruction.R"
    if destr_script.exists() and raw_outdir is not None:
        # Compare raw DA (native reference) against method DA
        for domain in ["cptac", "ccle"]:
            native_da = raw_outdir / "representation_da" / domain / "da_limma_result.csv"
            method_da = da_dir / domain / "da_limma_result.csv"
            if native_da.exists() and method_da.exists():
                logger.info("    Biology destruction: %s vs raw (%s)", method_key, domain)
                try:
                    run_r_script(destr_script, args=[
                        "--native-da", str(native_da),
                        "--method-da", str(method_da),
                        "--method", method_key,
                        "--outdir", str(calib_dir),
                    ], cwd=paths.root)
                except Exception as e:
                    logger.error("    Biology destruction failed (%s): %s", domain, e)

    # ── Marker sanity (per domain) ──────────────────────────────────────
    sanity_script = calibration_r_dir / "marker_sanity.R"
    if sanity_script.exists() and marker_info:
        for domain in ["cptac", "ccle"]:
            da_file = da_dir / domain / "da_limma_result.csv"
            if da_file.exists():
                logger.info("    Marker sanity: %s / %s", method_key, domain)
                try:
                    run_r_script(sanity_script, args=[
                        "--da-result", str(da_file),
                        "--markers", marker_info[0],
                        "--expected-signs", marker_info[1],
                        "--method", method_key,
                        "--domain", domain.upper(),
                        "--task", task.name,
                        "--outdir", str(calib_dir),
                    ], cwd=paths.root)
                except Exception as e:
                    logger.error("    Marker sanity failed (%s/%s): %s", method_key, domain, e)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Run benchmark evaluations")
    parser.add_argument("--config", default="configs/benchmark/default.yaml")
    parser.add_argument("--methods", nargs="+", default=None)
    parser.add_argument("--tasks", nargs="+", default=None)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--skip-calibration", action="store_true",
                        help="Skip calibration modules (permutation, ceiling, etc.)")
    parser.add_argument("--skip-validation", action="store_true",
                        help="Skip contrast direction validation")
    parser.add_argument(
        "--processed-dir",
        default="data/processed",
        help="Directory with shared_gene_matrix_<task>.csv and sample_meta_<task>.csv (from run_preprocessing.py)",
    )
    parser.add_argument(
        "--benchmark-results-dir",
        default="reports/benchmark_master/benchmark_results",
        help="Root directory for per-method benchmark outputs (avoids overwriting when comparing intersection vs union)",
    )
    args = parser.parse_args()

    bench_cfg = load_config(args.config)
    paths = ProjectPaths(args.repo_root)
    processed_dir = paths.resolve(args.processed_dir)
    benchmark_results_root = paths.resolve(args.benchmark_results_dir)

    task_names = args.tasks or bench_cfg.get("tasks", [])
    method_keys = args.methods or _expand_methods(bench_cfg.get("methods", []))

    logger.info("=" * 60)
    logger.info("  BENCHMARK RUNNER (with calibration)")
    logger.info("=" * 60)
    logger.info("Tasks: %s", ", ".join(task_names))
    logger.info("Methods: %s", ", ".join(method_keys))
    logger.info("Processed data: %s", processed_dir)
    logger.info("Benchmark results: %s", benchmark_results_root)
    logger.info("Calibration: %s", "SKIP" if args.skip_calibration else "ON")

    all_results = []
    import pandas as pd

    for task_name in task_names:
        task = TaskDefinition.from_yaml(task_name)
        logger.info("\n" + "=" * 50)
        logger.info("  Task: %s (%s)", task.name, task.contrast)
        logger.info("=" * 50)

        # Load metadata
        processed_meta = processed_dir / f"sample_meta_{task_name}.csv"
        if processed_meta.exists():
            sample_meta = pd.read_csv(processed_meta)
            logger.info("  Loaded preprocessed metadata: %d samples", len(sample_meta))
        else:
            from harmonize.preprocessing.subsets import build_subtype_subset, build_bvl_subset
            if task_name == "breast_subtype":
                sample_meta = build_subtype_subset(paths, task.raw_config)
            else:
                sample_meta = build_bvl_subset(paths, task.raw_config)
            sample_meta = sample_meta.drop_duplicates(subset=["sample_id"], keep="first")

        _, align_gene_index = _load_processed_matrix_if_present(
            paths, processed_dir, task.name, sample_meta,
        )

        # ── Phase A: Contrast validation ─────────────────────────────────
        if not args.skip_validation:
            validation_ok = _run_contrast_validation(task, paths, processed_dir)
            if not validation_ok:
                logger.error("  Contrast validation FAILED for %s — continuing with caution", task_name)

        # Load raw reference matrix for fixed-basis PCA
        raw_ref_matrix = _load_raw_reference(task, paths, sample_meta, processed_dir)

        raw_outdir = None  # Track raw method outdir for biology destruction comparison

        for method_key in method_keys:
            logger.info("\n  Method: %s", method_key)

            align = align_gene_index if method_key != "raw" else None
            result = _load_method_result(
                method_key, task, paths, sample_meta, processed_dir,
                align_gene_index=align,
            )
            if result is None:
                logger.warning("  Skipping %s — could not load", method_key)
                continue

            logger.info("  Loaded: %d genes x %d samples", result.n_genes, result.n_samples)

            outdir = benchmark_results_root / method_key / task_name
            outdir.mkdir(parents=True, exist_ok=True)

            if method_key == "raw":
                raw_outdir = outdir

            # ── Core benchmark (DA, markers, structure, matching) ─────────
            bench_result = run_benchmark(
                result, task, bench_cfg, outdir,
                reference_matrix=raw_ref_matrix,
            )
            all_results.append(bench_result)

            # ── Phase B: Per-method calibration ───────────────────────────
            if not args.skip_calibration:
                logger.info("  Phase B: Calibration for %s x %s", method_key, task_name)
                try:
                    _run_calibration_for_method(
                        method_key, task, paths, outdir, sample_meta, processed_dir,
                        raw_outdir=raw_outdir,
                    )
                except Exception as e:
                    logger.error("  Calibration failed for %s: %s", method_key, e)

    # ── Polished marker profile plots (R) ────────────────────────────────
    levels = bench_cfg.get("levels", {})
    if levels.get("marker_profiles", True):
        logger.info("\n" + "=" * 50)
        logger.info("  Generating polished marker profile plots (R)")
        logger.info("=" * 50)
        try:
            profile_script = paths.r_benchmark_dir / "run_polished_profile_plots.R"
            if profile_script.exists():
                run_r_script(profile_script, cwd=paths.root)
                logger.info("  Profile plots saved to reports/benchmark_master/marker_profiles/")
            else:
                logger.warning("  run_polished_profile_plots.R not found at %s", profile_script)
        except Exception as e:
            logger.error("  Marker profile plots failed: %s", e)

    # ── Phase C: Cross-method assembly ───────────────────────────────────
    if all_results:
        logger.info("\n" + "=" * 50)
        logger.info("  Phase C: Cross-method assembly")
        logger.info("=" * 50)

        summary_rows = []
        for r in all_results:
            row = {"method": r["method"], "task": r["task"]}
            if "representation_da" in r and isinstance(r["representation_da"], dict):
                row.update({f"da_{k}": v for k, v in r["representation_da"].items()
                            if not isinstance(v, (dict, list))})
            if "structure" in r and isinstance(r["structure"], dict):
                row.update({f"struct_{k}": v for k, v in r["structure"].items()
                            if not isinstance(v, (dict, list))})
            if "markers" in r and isinstance(r["markers"], dict):
                row["markers_present"] = r["markers"].get("n_present", 0)
            summary_rows.append(row)

        summary_df = pd.DataFrame(summary_rows)

        # Merge calibration results from disk
        from harmonize.reporting.summary_tables import build_method_comparison_table
        full_summary = build_method_comparison_table(benchmark_results_root)

        if not full_summary.empty:
            summary_df = full_summary

        benchmark_results_root.mkdir(parents=True, exist_ok=True)
        summary_path = benchmark_results_root / "comparison_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        logger.info("  Comparison summary saved: %s", summary_path)

        # Tiered summary
        try:
            from harmonize.reporting.summary_tables import annotate_metric_tiers
            tiered = annotate_metric_tiers(summary_df)
            tiered_path = benchmark_results_root / "comparison_summary_tiered.csv"
            tiered.to_csv(tiered_path, index=False)
            logger.info("  Tiered summary saved: %s", tiered_path)
        except Exception as e:
            logger.warning("  Tiered summary generation failed: %s", e)

    # ── Calibration figures (R + ggplot2) ────────────────────────────────
    if not args.skip_calibration:
        logger.info("\n" + "=" * 50)
        logger.info("  Generating calibration figures")
        logger.info("=" * 50)
        calib_fig_script = paths.r_benchmark_dir / "calibration_figures.R"
        if calib_fig_script.exists():
            try:
                run_r_script(calib_fig_script, args=["--repo-root", str(paths.root)],
                             cwd=paths.root)
                logger.info("  Calibration figures saved")
            except Exception as e:
                logger.error("  Calibration figures failed: %s", e)

    # ── Meeting slide figures ────────────────────────────────────────────
    meeting_fig_script = paths.r_benchmark_dir / "generate_meeting_figures.R"
    if meeting_fig_script.exists():
        logger.info("\n" + "=" * 50)
        logger.info("  Generating meeting slide figures")
        logger.info("=" * 50)
        try:
            run_r_script(meeting_fig_script, args=["--repo-root", str(paths.root)],
                         cwd=paths.root)
            logger.info("  Meeting figures saved")
        except Exception as e:
            logger.error("  Meeting figures failed: %s", e)

    logger.info("\n" + "=" * 60)
    logger.info("  BENCHMARK COMPLETE")
    logger.info("=" * 60)


def _expand_methods(method_names: list[str]) -> list[str]:
    """Expand method names (bridge_aware -> bridge_shift, bridge_scale)."""
    expanded = []
    for m in method_names:
        if m == "bridge_aware":
            expanded.extend(["bridge_shift", "bridge_scale"])
        else:
            expanded.append(m)
    return expanded


if __name__ == "__main__":
    main()
