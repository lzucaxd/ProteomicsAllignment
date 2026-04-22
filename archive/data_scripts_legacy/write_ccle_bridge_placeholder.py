#!/usr/bin/env python3
"""Placeholder PDF: CCLE bridge histogram not generated (missing qc_bridge outputs)."""
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "reports" / "benchmark_v1" / "diagnostics" / "ccle_bridge_histogram.pdf"


def main():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return 1
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        "CCLE bridge histogram — not generated in Benchmark v1 batch.\n\n"
        "No qc_bridge/*.tsv found under data/results/CCLE_corrected/.\n"
        "Run bridge QC on CCLE msstats_input.tsv to populate.",
        ha="center",
        va="center",
        fontsize=11,
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    plt.close()
    print("Wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
