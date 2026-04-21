#!/usr/bin/env python3
"""
Step 2: Copy union matrices to methods/raw, run bridge_shift / bridge_scale via R,
and Celligner on union splits. Outputs:
  data/processed/methods/{raw,bridge_shift,bridge_scale,celligner}/transformed_{task}.csv
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TASKS = {
    "breast_subtype": {
        "matrix": "shared_gene_matrix_breast_subtype.csv",
        "meta": "sample_meta_breast_subtype.csv",
    },
    "breast_vs_lung": {
        "matrix": "shared_gene_matrix_breast_vs_lung.csv",
        "meta": "sample_meta_breast_vs_lung.csv",
    },
}


def run_cmd(cmd: list[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    env = dict(**os.environ)
    # Some developer machines source renv/activate.R from ~/.Rprofile even when this
    # repo does not ship renv. That breaks plain `Rscript` subprocess calls.
    if not (REPO / "renv" / "activate.R").is_file():
        env.setdefault("R_PROFILE_USER", "/dev/null")
    subprocess.run(cmd, check=True, cwd=str(REPO), env=env)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=REPO)
    ap.add_argument(
        "--processed-dir",
        type=Path,
        default=REPO / "data" / "processed" / "union",
        help="Union matrices + sample_meta_*.csv",
    )
    ap.add_argument(
        "--out-methods-root",
        type=Path,
        default=REPO / "data" / "processed" / "methods",
    )
    ap.add_argument("--skip-celligner", action="store_true")
    args = ap.parse_args()

    repo = args.repo.resolve()
    proc = args.processed_dir.resolve()
    alt_union = repo / "data" / "processed" / "union"
    if not (proc / "shared_gene_matrix_breast_subtype.csv").is_file() and alt_union.is_dir():
        proc = alt_union
        print("Using union matrices from:", proc)
    root = args.out_methods_root.resolve()
    bridge_r = repo / "scripts" / "benchmark" / "bridge_aware_correction.R"
    celligner_py = repo / "scripts" / "benchmark" / "celligner_union_task.py"

    for m in ("raw", "bridge_shift", "bridge_scale", "celligner"):
        (root / m).mkdir(parents=True, exist_ok=True)

    for task, names in TASKS.items():
        mat = proc / names["matrix"]
        meta = proc / names["meta"]
        if not mat.is_file() or not meta.is_file():
            print(f"SKIP {task}: missing {mat} or {meta}", file=sys.stderr)
            continue

        raw_out = root / "raw" / f"transformed_{task}.csv"
        shutil.copy2(mat, raw_out)
        print(f"Copied raw -> {raw_out}")

        shift_out = root / "bridge_shift" / f"transformed_{task}.csv"
        run_cmd(
            [
                "Rscript",
                "--vanilla",
                str(bridge_r),
                "--union-task-matrix",
                "--repo",
                str(repo),
                "--matrix",
                str(mat),
                "--meta",
                str(meta),
                "--mode",
                "shift_only",
                "--out",
                str(shift_out),
            ]
        )

        scale_out = root / "bridge_scale" / f"transformed_{task}.csv"
        run_cmd(
            [
                "Rscript",
                "--vanilla",
                str(bridge_r),
                "--union-task-matrix",
                "--repo",
                str(repo),
                "--matrix",
                str(mat),
                "--meta",
                str(meta),
                "--mode",
                "shift_and_scale",
                "--out",
                str(scale_out),
            ]
        )

        if not args.skip_celligner:
            cg_out = root / "celligner" / f"transformed_{task}.csv"
            run_cmd(
                [
                    sys.executable,
                    str(celligner_py),
                    "--matrix",
                    str(mat),
                    "--meta",
                    str(meta),
                    "--out",
                    str(cg_out),
                    "--min-obs-frac",
                    "0",
                    "--no-umap",
                ]
            )

    print("regenerate_methods_union: done")


if __name__ == "__main__":
    main()
