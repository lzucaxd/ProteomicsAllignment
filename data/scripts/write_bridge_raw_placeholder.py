#!/usr/bin/env python3
"""Write a one-page PDF explaining raw-scale bridge plots were not generated in this batch."""
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "reports" / "benchmark_v1" / "diagnostics" / "bridge_boxplots_raw.pdf"

def main():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return 1
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")
    ax.text(
        0.5,
        0.55,
        "Raw-scale bridge channel diagnostics (not generated in Benchmark v1 batch)\n\n"
        "Per-(Run, Channel) linear intensities live in large msstats_input.tsv files.\n"
        "Use data/scripts/bridge_qc_cptac_ccle_same_scale.R or DuckDB aggregation\n"
        "to produce raw medians; current pack uses log2 summaries from qc_bridge/*.tsv.",
        ha="center",
        va="center",
        fontsize=11,
        wrap=True,
    )
    fig.savefig(OUT, dpi=150)
    plt.close()
    print("Wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
