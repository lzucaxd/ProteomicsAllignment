#!/usr/bin/env python3
"""Scatter CPTAC vs CCLE gene-level log2FC from shared_feature_table.csv."""

import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TABLE = REPO / "reports" / "benchmark_v1" / "shared_feature_table.csv"
OUT_DIR = REPO / "reports" / "benchmark_v1" / "diagnostics"


def main() -> int:
    if not TABLE.exists():
        print("Run build_benchmark_v1_artifacts.py first", file=sys.stderr)
        return 1
    xs, ys, mk = [], [], []
    markers = {"FOXA1", "GATA3", "EGFR", "KRT5", "KRT17", "ESR1", "FOXC1"}
    with TABLE.open(encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                xa = float(row["CPTAC_log2FC"])
                xb = float(row["CCLE_log2FC"])
            except (KeyError, ValueError):
                continue
            xs.append(xa)
            ys.append(xb)
            mk.append(row["GeneSymbol"] in markers)

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib required for PDF/PNG", file=sys.stderr)
        return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.axhline(0, color="gray", lw=0.6)
    ax.axvline(0, color="gray", lw=0.6)
    ax.scatter(
        [xs[i] for i in range(len(xs)) if not mk[i]],
        [ys[i] for i in range(len(ys)) if not mk[i]],
        s=8,
        alpha=0.35,
        c="steelblue",
        label="Other genes",
    )
    ax.scatter(
        [xs[i] for i in range(len(xs)) if mk[i]],
        [ys[i] for i in range(len(ys)) if mk[i]],
        s=36,
        alpha=0.9,
        c="darkorange",
        edgecolors="k",
        linewidths=0.4,
        label="Panel markers (subset)",
    )
    ax.set_xlabel("CPTAC log2FC (Luminal − Basal), gene-level collapsed")
    ax.set_ylabel("CCLE log2FC (Luminal − Basal), gene-level collapsed")
    ax.set_title("Shared-gene raw benchmark: cross-domain fold changes")
    ax.text(
        0.02,
        0.98,
        "Gene-level MSstatsTMT DA; one protein row per gene (min adj.P)",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        color="dimgray",
    )
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    base = OUT_DIR / "raw_fc_scatter"
    fig.savefig(base.with_suffix(".pdf"), dpi=200)
    fig.savefig(base.with_suffix(".png"), dpi=200)
    plt.close()
    print("Wrote", base.with_suffix(".pdf"), base.with_suffix(".png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
