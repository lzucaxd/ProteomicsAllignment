"""Celligner-aligned representation method wrapper.

Calls the existing Python script (run_celligner_all_data.py) which handles
Celligner fitting, alignment, PCA/UMAP generation.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from harmonize.methods.base import MethodInterface, MethodResult
from harmonize.utils.r_bridge import run_python_script
from harmonize.utils.paths import ProjectPaths
from harmonize.utils.io import load_gene_matrix

logger = logging.getLogger(__name__)


class CellignerMethod(MethodInterface):
    name = "celligner"
    display_name = "Celligner"

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

        python_cmd = config.get("python_cmd", "python3")
        conda_env = config.get("conda_env")
        if conda_env:
            python_cmd = f"/opt/anaconda3/envs/{conda_env}/bin/python"

        script = self.paths.resolve(
            config.get("python_script", "scripts/benchmark/run_celligner_all_data.py")
        )

        logger.info("Running Celligner alignment...")
        run_python_script(
            script,
            python_cmd=python_cmd,
            cwd=self.paths.root,
            timeout=7200,
        )

        return self.load_existing(sample_meta)

    def load_existing(self, sample_meta: pd.DataFrame | None = None) -> MethodResult:
        """Load pre-computed Celligner results."""
        cell_dir = self.paths.celligner_dir
        mat_path = cell_dir / "celligner_aligned_matrix.csv"
        meta_path = cell_dir / "sample_metadata.csv"

        if not mat_path.exists():
            raise FileNotFoundError(f"Celligner matrix not found: {mat_path}")

        # Celligner matrix is samples x genes, needs transposing
        raw = pd.read_csv(mat_path, index_col=0)
        matrix = raw.T
        matrix.index = matrix.index.astype(str)
        matrix.index.name = "GeneSymbol"

        if sample_meta is None and meta_path.exists():
            sample_meta = pd.read_csv(meta_path)
        elif sample_meta is None:
            sample_meta = pd.DataFrame()

        feature_meta = pd.DataFrame({
            "gene": matrix.index.tolist(),
            "included": True,
        })

        return MethodResult(
            matrix=matrix,
            sample_meta=sample_meta,
            feature_meta=feature_meta,
            method_name=self.name,
            display_name=self.display_name,
            notes="Celligner-aligned (cPCA + MNN). Values are z-scored.",
            transforms_values=True,
            value_scale="z_scored",
            qc_paths=[
                str(cell_dir / "pca_pre.png"),
                str(cell_dir / "pca_post.png"),
                str(cell_dir / "umap_post.png"),
            ],
        )
