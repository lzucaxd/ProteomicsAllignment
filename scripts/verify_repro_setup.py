#!/usr/bin/env python3
"""
Verify a clean clone is ready to run preprocessing + benchmark Python/R steps.

Usage (from repo root):
  python3 scripts/verify_repro_setup.py
  python3 scripts/verify_repro_setup.py --require-data
  python3 scripts/verify_repro_setup.py --skip-r   # Python stack only (e.g. CI)

--require-data checks paths in configs/preprocessing/default.yaml (gene matrices).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)


def _ok(msg: str) -> None:
    print(f"OK   {msg}")


def check_python_version() -> bool:
    if sys.version_info < (3, 10):
        _fail(f"Python >= 3.10 required (harmonize); got {sys.version_info.major}.{sys.version_info.minor}")
        return False
    _ok(f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
    return True


def check_imports(repo: Path) -> bool:
    sys.path.insert(0, str(repo / "src"))
    try:
        import yaml  # noqa: F401
        import numpy  # noqa: F401
        import pandas  # noqa: F401
        import scipy  # noqa: F401
        import sklearn  # noqa: F401
        import statsmodels  # noqa: F401
        from harmonize.utils.config import load_config  # noqa: F401

        _ok("Python imports (harmonize + scientific stack)")
        return True
    except ImportError as e:
        _fail(f"Import error: {e}. Install deps: pip install -r requirements.txt or pip install -e .")
        return False


def check_rscript() -> bool:
    try:
        r = subprocess.run(
            ["Rscript", "--version"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if r.returncode != 0:
            _fail("Rscript returned non-zero")
            return False
        line = (r.stdout or r.stderr or "").strip().splitlines()
        ver = line[0] if line else "(version unknown)"
        _ok(f"Rscript available ({ver})")
        return True
    except FileNotFoundError:
        _fail("Rscript not found on PATH (install R >= 4.1 and use `Rscript install_r_packages.R`)")
        return False
    except subprocess.TimeoutExpired:
        _fail("Rscript --version timed out")
        return False


def check_data_paths(repo: Path) -> bool:
    try:
        sys.path.insert(0, str(repo / "src"))
        from harmonize.utils.config import load_config
    except ImportError:
        _fail("Cannot load harmonize; fix imports before --require-data")
        return False

    cfg_path = repo / "configs" / "preprocessing" / "default.yaml"
    if not cfg_path.is_file():
        _fail(f"Missing {cfg_path}")
        return False

    cfg = load_config(cfg_path)
    ds = cfg.get("data_sources", {})
    ok = True
    for study, info in ds.get("cptac", {}).items():
        rel = info.get("gene_matrix")
        if not rel:
            continue
        p = repo / rel
        if p.is_file():
            _ok(f"CPTAC gene_matrix {study}: {rel}")
        else:
            _fail(f"Missing CPTAC gene_matrix {study}: {p}")
            ok = False

    ccle = ds.get("ccle", {})
    for key in ("gene_matrix", "sample_info"):
        rel = ccle.get(key)
        if not rel:
            continue
        p = repo / rel
        if p.is_file():
            _ok(f"CCLE {key}: {rel}")
        else:
            _fail(f"Missing CCLE {key}: {p}")
            ok = False

    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--require-data",
        action="store_true",
        help="Require gene_matrix.csv paths from configs/preprocessing/default.yaml",
    )
    parser.add_argument(
        "--skip-r",
        action="store_true",
        help="Do not require Rscript on PATH (full benchmark still needs R).",
    )
    args = parser.parse_args()

    repo = _repo_root()
    print(f"Repository root: {repo}\n", flush=True)

    checks = [
        check_python_version(),
        check_imports(repo),
    ]
    if not args.skip_r:
        checks.append(check_rscript())
    if args.require_data:
        print("")
        checks.append(check_data_paths(repo))

    if all(checks):
        print("\nAll checks passed.", flush=True)
        if args.skip_r:
            print("(R was skipped; run without --skip-r locally before overnight benchmark.)", flush=True)
        if not args.require_data:
            print("(Re-run with --require-data after placing CPTAC/CCLE gene matrices.)", flush=True)
        return 0

    print("\nOne or more checks failed.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
