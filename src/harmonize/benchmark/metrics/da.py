"""Representation-level differential abundance via limma (R subprocess)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from harmonize.utils.r_bridge import run_r_script
from harmonize.utils.paths import ProjectPaths

logger = logging.getLogger(__name__)


def run_limma_da(
    matrix: pd.DataFrame,
    sample_meta: pd.DataFrame,
    task_name: str,
    method_name: str,
    outdir: str | Path,
    paths: ProjectPaths | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Run representation-level limma DA by calling existing R evaluation helpers.

    Dispatches to scripts/benchmark/benchmark_runner.R which sources
    evaluation_helpers.R and the appropriate task module.

    Returns dict of per-domain DA result DataFrames.
    """
    if paths is None:
        paths = ProjectPaths()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Save inputs for R
    matrix.to_csv(outdir / "_input_matrix.csv")
    sample_meta.to_csv(outdir / "_input_meta.csv", index=False)

    runner = paths.r_benchmark_dir / "benchmark_runner.R"
    run_r_script(
        runner,
        args=[
            "--matrix", str(outdir / "_input_matrix.csv"),
            "--meta", str(outdir / "_input_meta.csv"),
            "--task", task_name,
            "--method", method_name,
            "--outdir", str(outdir),
        ],
        cwd=paths.root,
    )

    # Load results
    results = {}
    for domain_dir in outdir.iterdir():
        if domain_dir.is_dir():
            da_file = domain_dir / "da_limma_result.csv"
            if da_file.exists():
                results[domain_dir.name] = pd.read_csv(da_file)

    return results


def run_limma_da_python(
    matrix: pd.DataFrame,
    groups: pd.Series,
    contrast_name: str = "GroupB_vs_GroupA",
) -> pd.DataFrame:
    """
    DEPRECATED — Python-native two-group Welch t-tests.

    WARNING: This function determines group order from ``groups.unique()``
    which depends on metadata row order, NOT alphabetical sorting. This
    causes sign-flipped logFC when CPTAC and CCLE metadata list conditions
    in different orders. The benchmark pipeline now uses R limma via
    ``limma_da_wrapper.R`` with explicit contrast levels. This fallback is
    retained only for environments without R.
    """
    import warnings
    warnings.warn(
        "run_limma_da_python uses metadata-order-dependent contrast direction. "
        "Use R limma via runner._run_limma_r() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    from scipy import stats
    import numpy as np

    group_a = groups.unique()[0]
    group_b = groups.unique()[1] if len(groups.unique()) > 1 else group_a
    mask_a = groups == group_a
    mask_b = groups == group_b
    cols_a = groups.index[mask_a]
    cols_b = groups.index[mask_b]

    results = []
    for gene in matrix.index:
        vals_a = matrix.loc[gene, cols_a].dropna().values
        vals_b = matrix.loc[gene, cols_b].dropna().values
        if len(vals_a) < 2 or len(vals_b) < 2:
            continue
        t_stat, p_val = stats.ttest_ind(vals_b, vals_a, equal_var=False)
        log_fc = np.nanmean(vals_b) - np.nanmean(vals_a)
        results.append({
            "gene": gene,
            "logFC": log_fc,
            "t": t_stat,
            "P.Value": p_val,
            "contrast": contrast_name,
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        from statsmodels.stats.multitest import multipletests
        df["adj.P.Val"] = multipletests(df["P.Value"], method="fdr_bh")[1]
    return df
