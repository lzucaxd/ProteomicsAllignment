"""Main benchmark orchestrator: runs all evaluation levels for a method x task."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from harmonize.benchmark.tasks import TaskDefinition
from harmonize.methods.base import MethodResult

logger = logging.getLogger(__name__)


def run_benchmark(
    method_result: MethodResult,
    task: TaskDefinition,
    benchmark_config: dict[str, Any],
    outdir: str | Path,
    reference_matrix: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """
    Run all enabled benchmark levels for one method on one task.

    Parameters
    ----------
    reference_matrix : optional raw/unharmonized matrix for fixed-basis PCA.

    Returns a summary dict with per-level results.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    results = {"method": method_result.method_name, "task": task.name}

    levels = benchmark_config.get("levels", {})

    # ── Level 2: Representation-level DA ────────────────────────────
    if levels.get("representation_level_da", True):
        logger.info("  Level 2: Representation-level DA")
        da_dir = outdir / "representation_da"
        try:
            da_results = _run_representation_da(method_result, task, da_dir)
            results["representation_da"] = da_results
        except Exception as e:
            logger.error("  Level 2 failed: %s", e)
            results["representation_da"] = {"error": str(e)}

    # ── Level 3: Marker profiles ────────────────────────────────────
    if levels.get("marker_profiles", True):
        logger.info("  Level 3: Marker profiles")
        try:
            marker_info = _check_markers(method_result, task)
            results["markers"] = marker_info
        except Exception as e:
            logger.error("  Level 3 failed: %s", e)

    # ── Level 4: Structure metrics + plots ──────────────────────────
    if levels.get("structure_metrics", True):
        logger.info("  Level 4: Structure metrics")
        struct_dir = outdir / "structure"
        struct_dir.mkdir(parents=True, exist_ok=True)
        try:
            from harmonize.benchmark.metrics.structure import compute_structure_metrics
            struct = compute_structure_metrics(
                method_result.matrix,
                method_result.sample_meta,
                config=benchmark_config.get("structure_metrics", {}),
                reference_matrix=reference_matrix,
            )
            results["structure"] = struct

            struct_df = pd.DataFrame([struct])
            struct_df.to_csv(struct_dir / "structure_summary.csv", index=False)
        except Exception as e:
            logger.error("  Level 4 metrics failed: %s", e)
            results["structure"] = {"error": str(e)}

        # PCA / UMAP plots
        try:
            from harmonize.benchmark.plots.structure import plot_pca_structure, plot_umap_structure
            plot_pca_structure(
                method_result.matrix, method_result.sample_meta,
                method_result.method_name, task.name, struct_dir,
            )
            sc = benchmark_config.get("structure_metrics", {})
            if sc.get("compute_umap", False):
                plot_umap_structure(
                    method_result.matrix, method_result.sample_meta,
                    method_result.method_name, task.name, struct_dir,
                    n_neighbors=sc.get("umap_n_neighbors", 15),
                    min_dist=sc.get("umap_min_dist", 0.3),
                )
        except Exception as e:
            logger.error("  Level 4 plots failed: %s", e)

    # ── Level 5: Matching metrics ───────────────────────────────────
    if levels.get("matching_metrics", False):
        logger.info("  Level 5: Matching metrics")
        match_dir = outdir / "matching"
        match_dir.mkdir(parents=True, exist_ok=True)
        try:
            from harmonize.benchmark.metrics.matching import compute_matching_metrics
            match = compute_matching_metrics(
                method_result.matrix,
                method_result.sample_meta,
                config=benchmark_config.get("matching_metrics", {}),
            )
            results["matching"] = match
        except Exception as e:
            logger.error("  Level 5 failed: %s", e)
            results["matching"] = {"error": str(e)}

    # ── Save summary ────────────────────────────────────────────────
    _save_summary(results, outdir)
    return results


def _run_representation_da(
    method_result: MethodResult,
    task: TaskDefinition,
    outdir: Path,
) -> dict:
    """Run per-domain limma DA (via R subprocess) and cross-domain agreement."""
    from harmonize.benchmark.metrics.agreement import compute_fc_agreement, summarize_agreement

    outdir.mkdir(parents=True, exist_ok=True)
    meta = method_result.sample_meta
    mat = method_result.matrix

    if meta.empty or "domain" not in meta.columns:
        return {"error": "No sample metadata with domain info"}

    valid_cols = [c for c in mat.columns if c in set(meta["sample_id"])]
    if len(valid_cols) == 0:
        return {"error": "No matching samples between matrix and metadata"}

    # Parse contrast into explicit levels: "Luminal_vs_Basal" -> ("Basal", "Luminal")
    contrast_parts = task.contrast.split("_vs_")
    if len(contrast_parts) == 2:
        contrast_b, contrast_a = contrast_parts[0], contrast_parts[1]
    else:
        contrast_a, contrast_b = sorted(meta["condition"].dropna().unique())[:2]

    da_per_domain = _run_limma_r(mat, meta, valid_cols, contrast_a, contrast_b,
                                  task.contrast, outdir)

    summary = {"n_domains_with_da": len(da_per_domain)}
    if "CPTAC" in da_per_domain and "CCLE" in da_per_domain:
        agreement = compute_fc_agreement(da_per_domain["CPTAC"], da_per_domain["CCLE"])
        agreement.to_csv(outdir / "fc_agreement.csv", index=False)
        summary.update(summarize_agreement(agreement))

        try:
            from harmonize.benchmark.plots.diagnostics import plot_fc_scatter
            plot_fc_scatter(
                da_per_domain["CPTAC"], da_per_domain["CCLE"],
                method_result.method_name, task.name, outdir,
                markers=task.markers,
            )
        except Exception as e:
            logger.warning("  FC scatter plot failed: %s", e)

    return summary


def _run_limma_r(
    mat: pd.DataFrame,
    meta: pd.DataFrame,
    valid_cols: list[str],
    contrast_a: str,
    contrast_b: str,
    contrast_name: str,
    outdir: Path,
) -> dict[str, pd.DataFrame]:
    """Call the R limma wrapper and read back per-domain results."""
    from harmonize.utils.r_bridge import run_r_script
    from harmonize.utils.paths import ProjectPaths

    paths = ProjectPaths()

    tmp_matrix = outdir / "_input_matrix.csv"
    tmp_meta = outdir / "_input_meta.csv"
    mat[valid_cols].to_csv(tmp_matrix)
    meta[meta["sample_id"].isin(valid_cols)].to_csv(tmp_meta, index=False)

    wrapper_script = (
        paths.root / "src" / "harmonize" / "benchmark" / "calibration" / "limma_da_wrapper.R"
    )

    run_r_script(
        wrapper_script,
        args=[
            "--matrix", str(tmp_matrix),
            "--meta", str(tmp_meta),
            "--contrast-a", contrast_a,
            "--contrast-b", contrast_b,
            "--contrast-name", contrast_name,
            "--outdir", str(outdir),
        ],
        cwd=paths.root,
    )

    da_per_domain = {}
    for domain in ["CPTAC", "CCLE"]:
        da_file = outdir / domain.lower() / "da_limma_result.csv"
        if da_file.exists():
            da_per_domain[domain] = pd.read_csv(da_file)
            logger.info("  %s limma DA: %d genes", domain, len(da_per_domain[domain]))

    # Clean up temp files
    for f in [tmp_matrix, tmp_meta]:
        try:
            f.unlink()
        except OSError:
            pass

    return da_per_domain


def _check_markers(method_result: MethodResult, task: TaskDefinition) -> dict:
    """Check marker availability and directions."""
    from harmonize.benchmark.metrics.markers import check_marker_availability

    avail = check_marker_availability(method_result.matrix, task.markers)
    n_present = sum(avail.values())
    return {
        "n_requested": len(task.markers),
        "n_present": n_present,
        "present": [m for m, v in avail.items() if v],
        "absent": [m for m, v in avail.items() if not v],
    }


def _save_summary(results: dict, outdir: Path) -> None:
    """Save benchmark summary to a text file."""
    outdir.mkdir(parents=True, exist_ok=True)
    lines = [f"Benchmark: {results['method']} x {results['task']}", ""]
    for key, val in results.items():
        if key in ("method", "task"):
            continue
        lines.append(f"[{key}]")
        if isinstance(val, dict):
            for k, v in val.items():
                lines.append(f"  {k}: {v}")
        else:
            lines.append(f"  {val}")
        lines.append("")
    (outdir / "benchmark_summary.txt").write_text("\n".join(lines))
