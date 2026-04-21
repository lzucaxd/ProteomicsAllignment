"""Bridge-aware per-protein harmonization method wrapper.

Calls existing R scripts:
  1. extract_bridge_summaries.R  — extract bridge channel data
  2. bridge_aware_correction.R   — apply shift / shift+scale correction
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from harmonize.methods.base import MethodInterface, MethodResult
from harmonize.utils.r_bridge import run_r_script
from harmonize.utils.paths import ProjectPaths
from harmonize.utils.io import load_gene_matrix

logger = logging.getLogger(__name__)


class BridgeAwareMethod(MethodInterface):
    name = "bridge_aware"
    display_name = "Bridge-Aware"

    def __init__(self, paths: ProjectPaths | None = None):
        self.paths = paths or ProjectPaths()

    def run(
        self,
        cptac_matrix: pd.DataFrame,
        ccle_matrix: pd.DataFrame,
        sample_meta: pd.DataFrame,
        config: dict[str, Any],
        outdir: str | Path,
    ) -> MethodResult:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Step 1: extract bridge summaries (if not already present)
        bridge_dir = self.paths.bridge_dir
        bridge_cptac = bridge_dir / "bridge_summary_cptac.tsv"
        bridge_ccle = bridge_dir / "bridge_summary_ccle.tsv"

        if not bridge_cptac.exists() or not bridge_ccle.exists():
            logger.info("Extracting bridge summaries...")
            extract_script = self.paths.r_benchmark_dir / "extract_bridge_summaries.R"
            run_r_script(extract_script, cwd=self.paths.root)

        # Step 2: run bridge correction
        logger.info("Running bridge-aware correction...")
        correction_script = self.paths.r_benchmark_dir / "bridge_aware_correction.R"
        run_r_script(correction_script, cwd=self.paths.root)

        # Load outputs for each mode
        modes = config.get("modes", ["shift_only", "shift_and_scale"])
        results = {}

        for mode in modes:
            if mode == "shift_only":
                mat_path = bridge_dir / "bridge_aware_shift_only_matrix.csv"
                suffix = "shift"
            else:
                mat_path = bridge_dir / "bridge_aware_shift_scale_matrix.csv"
                suffix = "scale"

            if not mat_path.exists():
                logger.warning("Bridge %s matrix not found: %s", mode, mat_path)
                continue

            matrix = load_gene_matrix(mat_path)
            valid_sids = set(sample_meta["sample_id"])
            matrix = matrix[[c for c in matrix.columns if c in valid_sids]]

            feature_meta = pd.DataFrame({
                "gene": matrix.index.tolist(),
                "included": True,
                "bridge_corrected": True,
            })

            result = MethodResult(
                matrix=matrix,
                sample_meta=sample_meta,
                feature_meta=feature_meta,
                method_name=f"bridge_{suffix}",
                display_name=f"Bridge {'Shift-Only' if mode == 'shift_only' else 'Shift+Scale'}",
                notes=f"Bridge-aware {mode} correction. See {bridge_dir}/bridge_aware_{suffix}*_qc.md",
                transforms_values=True,
                value_scale="log2_abundance_calibrated",
                qc_paths=[str(bridge_dir / f"bridge_aware_{suffix}_qc.md")],
            )
            results[mode] = result

        # Return shift_only as primary; save both
        primary = results.get("shift_only") or results.get("shift_and_scale")
        if primary is None:
            raise RuntimeError("No bridge-aware matrices produced")

        primary.save(outdir)
        return primary

    def load_existing(self, mode: str = "shift_only") -> MethodResult | None:
        """Load pre-computed bridge results without re-running."""
        bridge_dir = self.paths.bridge_dir
        if mode == "shift_only":
            mat_path = bridge_dir / "bridge_aware_shift_only_matrix.csv"
            name, display = "bridge_shift", "Bridge Shift-Only"
        else:
            mat_path = bridge_dir / "bridge_aware_shift_scale_matrix.csv"
            name, display = "bridge_scale", "Bridge Shift+Scale"

        if not mat_path.exists():
            return None

        matrix = load_gene_matrix(mat_path)
        return MethodResult(
            matrix=matrix,
            sample_meta=pd.DataFrame(),
            feature_meta=pd.DataFrame({"gene": matrix.index.tolist(), "included": True}),
            method_name=name,
            display_name=display,
            transforms_values=True,
            value_scale="log2_abundance_calibrated",
        )
