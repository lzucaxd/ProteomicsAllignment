"""Centralized path management for the benchmark framework."""

from __future__ import annotations

from pathlib import Path
from typing import Any


class ProjectPaths:
    """Resolves all project paths from config and creates output directories on demand."""

    def __init__(self, repo_root: str | Path | None = None, config: dict[str, Any] | None = None):
        if repo_root is None:
            repo_root = Path(__file__).resolve().parents[3]
        self.root = Path(repo_root).resolve()
        self._cfg = config or {}

    def resolve(self, relative: str) -> Path:
        """Resolve a path relative to repo root."""
        return self.root / relative

    def ensure_dir(self, relative: str) -> Path:
        """Resolve and create a directory."""
        p = self.resolve(relative)
        p.mkdir(parents=True, exist_ok=True)
        return p

    # ── Data sources ────────────────────────────────────────────────
    @property
    def cptac_breast_matrix(self) -> Path:
        return self.resolve("data/results/PDC000120/gene_matrix.csv")

    @property
    def cptac_lung_matrix(self) -> Path:
        return self.resolve("data/results/PDC000153/gene_matrix.csv")

    @property
    def cptac_ovarian_matrix(self) -> Path:
        return self.resolve("data/results/PDC000127/gene_matrix.csv")

    @property
    def cptac_uterine_matrix(self) -> Path:
        return self.resolve("data/results/PDC000204/gene_matrix.csv")

    @property
    def ccle_matrix(self) -> Path:
        return self.resolve("data/results/CCLE_corrected/gene_matrix.csv")

    @property
    def ccle_sample_info(self) -> Path:
        return self.resolve("data/ccle_peptide/sample_info_ccle.csv")

    @property
    def subtype_mapping(self) -> Path:
        return self.resolve(
            "data/annotations/cptac/PDC000120/gene_matrix_subtype_mapping.csv"
        )

    # ── Processed outputs ───────────────────────────────────────────
    @property
    def processed_dir(self) -> Path:
        return self.ensure_dir("data/processed")

    # ── Method outputs ──────────────────────────────────────────────
    @property
    def bridge_dir(self) -> Path:
        return self.resolve("reports/benchmark_master/methods/bridge_aware")

    @property
    def celligner_dir(self) -> Path:
        return self.resolve("reports/benchmark_master/celligner_all")

    @property
    def raw_method_dir(self) -> Path:
        return self.resolve("reports/benchmark_master/methods/raw")

    # ── Benchmark outputs ───────────────────────────────────────────
    @property
    def benchmark_root(self) -> Path:
        return self.resolve("reports/benchmark_master")

    def benchmark_output(self, method: str, task: str, level: str) -> Path:
        return self.ensure_dir(f"reports/benchmark_master/{level}/{method}/{task}")

    # ── R scripts ───────────────────────────────────────────────────
    @property
    def r_methods_dir(self) -> Path:
        return self.resolve("scripts/methods")

    @property
    def r_benchmark_dir(self) -> Path:
        return self.resolve("scripts/benchmark")

    # ── Configs ─────────────────────────────────────────────────────
    @property
    def configs_dir(self) -> Path:
        return self.resolve("configs")
