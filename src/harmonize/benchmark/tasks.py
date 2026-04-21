"""Task definitions loaded from YAML configs."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from harmonize.utils.config import load_task_config


@dataclass
class TaskDefinition:
    """A benchmark task specification."""

    name: str
    contrast: str
    block_order: list[str]
    markers: list[str]
    expected_directions: dict[str, str] = field(default_factory=dict)
    feature_level: str = "gene"
    native_inference: dict[str, str] = field(default_factory=dict)
    representation_inference: str = "limma"
    raw_config: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, task_name: str, config_dir: str | None = None) -> TaskDefinition:
        """Load a TaskDefinition from a YAML config file."""
        cfg = load_task_config(task_name, config_dir)
        return cls(
            name=cfg["task_name"],
            contrast=cfg.get("contrast", ""),
            block_order=cfg.get("block_order", []),
            markers=cfg.get("markers", []),
            expected_directions=cfg.get("expected_marker_directions", {}),
            feature_level=cfg.get("feature_level", "gene"),
            native_inference=cfg.get("native_domain_inference", {}),
            representation_inference=cfg.get("representation_level_inference", "limma"),
            raw_config=cfg,
        )

    def validate_sample_meta(self, meta: pd.DataFrame) -> list[str]:
        """Check that sample metadata has required columns. Returns list of warnings."""
        warnings = []
        required = {"sample_id", "domain", "condition"}
        missing = required - set(meta.columns)
        if missing:
            warnings.append(f"Missing columns: {missing}")

        domains = set(meta["domain"].unique()) if "domain" in meta.columns else set()
        if "CPTAC" not in domains:
            warnings.append("No CPTAC samples in metadata")
        if "CCLE" not in domains:
            warnings.append("No CCLE samples in metadata")

        return warnings
