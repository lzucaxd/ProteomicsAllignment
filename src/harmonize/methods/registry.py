"""Method registry: discover and instantiate harmonization methods by name."""

from __future__ import annotations

from typing import Any

from harmonize.methods.base import MethodInterface, MethodResult
from harmonize.methods.raw import RawMethod
from harmonize.methods.bridge_aware import BridgeAwareMethod
from harmonize.methods.celligner import CellignerMethod
from harmonize.utils.paths import ProjectPaths

_REGISTRY: dict[str, type[MethodInterface]] = {
    "raw": RawMethod,
    "bridge_aware": BridgeAwareMethod,
    "celligner": CellignerMethod,
}


def register_method(name: str, cls: type[MethodInterface]) -> None:
    """Register a new method class (for future extensions)."""
    _REGISTRY[name] = cls


def list_methods() -> list[str]:
    """Return names of all registered methods."""
    return list(_REGISTRY.keys())


def get_method(name: str, paths: ProjectPaths | None = None) -> MethodInterface:
    """Instantiate a method by name."""
    if name not in _REGISTRY:
        raise KeyError(
            f"Unknown method '{name}'. Available: {list(_REGISTRY.keys())}"
        )
    cls = _REGISTRY[name]
    if paths and hasattr(cls, "__init__"):
        try:
            return cls(paths=paths)
        except TypeError:
            return cls()
    return cls()


def run_method(
    name: str,
    config: dict[str, Any],
    paths: ProjectPaths | None = None,
    **kwargs,
) -> MethodResult:
    """Instantiate and run a method by name."""
    method = get_method(name, paths)
    return method.run(**kwargs, config=config)
