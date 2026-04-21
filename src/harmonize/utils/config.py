"""YAML config loading with path resolution."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML config file and return as a dict."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with open(path) as f:
        return yaml.safe_load(f) or {}


def load_task_config(task_name: str, config_dir: str | Path | None = None) -> dict[str, Any]:
    """Load a task config by name from configs/tasks/."""
    if config_dir is None:
        config_dir = Path("configs/tasks")
    return load_config(Path(config_dir) / f"{task_name}.yaml")


def load_method_config(method_name: str, config_dir: str | Path | None = None) -> dict[str, Any]:
    """Load a method config by name from configs/methods/."""
    if config_dir is None:
        config_dir = Path("configs/methods")
    return load_config(Path(config_dir) / f"{method_name}.yaml")
