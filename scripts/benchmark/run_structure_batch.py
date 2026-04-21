#!/usr/bin/env python3
"""
Compute structure metrics (full gene matrix) for each method × task; fixed-basis PCA from raw.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from harmonize.benchmark.metrics.structure import compute_structure_metrics  # noqa: E402


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=REPO)
    ap.add_argument("--methods-root", type=Path, default=REPO / "data" / "processed" / "methods")
    ap.add_argument("--meta-dir", type=Path, default=REPO / "data" / "processed_union")
    ap.add_argument("--results-root", type=Path, default=REPO / "reports" / "benchmark_master" / "benchmark_results")
    args = ap.parse_args()
    repo = args.repo.resolve()
    methods_root = args.methods_root.resolve()
    meta_dir = args.meta_dir.resolve()
    results_root = args.results_root.resolve()

    tasks = {
        "breast_subtype": "sample_meta_breast_subtype.csv",
        "breast_vs_lung": "sample_meta_breast_vs_lung.csv",
    }
    methods = ["raw", "bridge_shift", "bridge_scale", "celligner"]
    struct_cfg = {"n_pcs": 20, "knn_k": 15, "compute_umap": False}

    for task, meta_name in tasks.items():
        raw_mat_path = methods_root / "raw" / f"transformed_{task}.csv"
        if not raw_mat_path.is_file():
            print(f"SKIP structure {task}: missing {raw_mat_path}", file=sys.stderr)
            continue
        ref = load_matrix(raw_mat_path)
        meta_path = meta_dir / meta_name
        if not meta_path.is_file():
            for alt in (
                repo / "data" / "processed" / "union" / meta_name,
                repo / "data" / "processed_union" / meta_name,
            ):
                if alt.is_file():
                    meta_path = alt
                    break
        meta = pd.read_csv(meta_path)

        for m in methods:
            mat_path = methods_root / m / f"transformed_{task}.csv"
            if not mat_path.is_file():
                continue
            mat = load_matrix(mat_path)
            outdir = results_root / m / task / "structure"
            outdir.mkdir(parents=True, exist_ok=True)
            struct = compute_structure_metrics(
                mat, meta, config=struct_cfg, reference_matrix=ref
            )
            pd.DataFrame([struct]).to_csv(outdir / "structure_summary.csv", index=False)
            print(f"structure {m} {task} -> {outdir / 'structure_summary.csv'}")


if __name__ == "__main__":
    main()
