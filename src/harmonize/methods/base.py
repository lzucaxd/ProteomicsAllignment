"""Base interface and result type for all harmonization methods."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class MethodResult:
    """Standard return type from any harmonization method.

    All methods must produce this structure so the benchmark layer
    can consume them uniformly.
    """
    matrix: pd.DataFrame           # genes x samples
    sample_meta: pd.DataFrame      # sample_id, domain, condition, ...
    feature_meta: pd.DataFrame     # gene, included, ...
    method_name: str
    display_name: str = ""
    notes: str = ""
    qc_paths: list[str] = field(default_factory=list)
    transforms_values: bool = False
    value_scale: str = "log2_abundance"  # or "z_scored", "residual", etc.

    def __post_init__(self):
        if not self.display_name:
            self.display_name = self.method_name

    @property
    def n_genes(self) -> int:
        return len(self.matrix)

    @property
    def n_samples(self) -> int:
        return len(self.matrix.columns)

    def save(self, outdir: str | Path) -> None:
        """Save all components to disk."""
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.matrix.to_csv(outdir / "transformed_matrix.csv")
        self.sample_meta.to_csv(outdir / "sample_metadata.csv", index=False)
        self.feature_meta.to_csv(outdir / "feature_metadata.csv", index=False)
        (outdir / "method_notes.txt").write_text(self.notes)

    @classmethod
    def load(cls, outdir: str | Path, method_name: str = "unknown") -> MethodResult:
        """Load a previously saved MethodResult."""
        outdir = Path(outdir)
        matrix = pd.read_csv(outdir / "transformed_matrix.csv", index_col=0)
        matrix.index = matrix.index.astype(str)
        sample_meta = pd.read_csv(outdir / "sample_metadata.csv")
        feature_meta = pd.read_csv(outdir / "feature_metadata.csv")
        notes = ""
        notes_path = outdir / "method_notes.txt"
        if notes_path.exists():
            notes = notes_path.read_text()
        return cls(
            matrix=matrix,
            sample_meta=sample_meta,
            feature_meta=feature_meta,
            method_name=method_name,
            notes=notes,
        )


class MethodInterface(ABC):
    """Abstract base for all harmonization method wrappers."""

    name: str = "abstract"
    display_name: str = "Abstract"

    @abstractmethod
    def run(
        self,
        cptac_matrix: pd.DataFrame,
        ccle_matrix: pd.DataFrame,
        sample_meta: pd.DataFrame,
        config: dict[str, Any],
        outdir: str | Path,
    ) -> MethodResult:
        """Run the method and return a standardized result."""
        ...

    def load_result(self, outdir: str | Path) -> MethodResult:
        """Load a previously generated result."""
        return MethodResult.load(outdir, self.name)
