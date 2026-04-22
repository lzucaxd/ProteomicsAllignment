#!/usr/bin/env python3
"""
Subsample msstats_input.tsv rows: bridge (Condition==Norm) vs sample (else).
Produces separate figures for raw linear intensity and log2(intensity+1).
CPTAC and CCLE side-by-side in one PDF each.

--max-rows caps rows read per file (default 800000) for reproducible runtime.

Outputs:
  reports/benchmark_v1/diagnostics_feedback/figures/bridge_boxplots_raw.pdf
  reports/benchmark_v1/diagnostics_feedback/figures/bridge_boxplots_log.pdf
  sample_boxplots_raw.pdf / sample_boxplots_log.pdf — copies of same figures (bridge vs sample is the comparison).
"""
from __future__ import annotations

import csv
import math
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "reports" / "benchmark_v1" / "diagnostics_feedback" / "figures"


def log2p1(x: float) -> float:
    return math.log2(max(x, 0.0) + 1.0)


def collect(path: Path, max_rows: int) -> tuple[list[float], list[float], list[float], list[float]]:
    raw_n, raw_s, log_n, log_s = [], [], [], []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(r):
            if i >= max_rows:
                break
            try:
                intensity = float(row["Intensity"])
            except (KeyError, ValueError):
                continue
            cond = (row.get("Condition") or "").strip()
            if cond == "Norm":
                raw_n.append(intensity)
                log_n.append(log2p1(intensity))
            else:
                raw_s.append(intensity)
                log_s.append(log2p1(intensity))
    return raw_n, raw_s, log_n, log_s


def plot_four_panels(
    cptac_raw_n,
    cptac_raw_s,
    ccle_raw_n,
    ccle_raw_s,
    ylab: str,
    title: str,
    stem: str,
    max_rows: int,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib required", file=sys.stderr)
        return

    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    def one(ax, a, b, lab):
        ax.boxplot([a, b], tick_labels=["Bridge/Norm", "Sample"], showfliers=False)
        ax.set_title(lab)
        ax.set_ylabel(ylab)

    one(axes[0], cptac_raw_n, cptac_raw_s, "CPTAC")
    one(axes[1], ccle_raw_n, ccle_raw_s, "CCLE")
    fig.suptitle(
        title + f"\n(First {max_rows:,} msstats_input rows per cohort; see raw_vs_log_boxplot_notes.md)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.pdf", dpi=200)
    fig.savefig(OUT / f"{stem}.png", dpi=200)
    plt.close()


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--max-rows", type=int, default=800_000)
    args = ap.parse_args()

    cptac_p = REPO / "data" / "results" / "PDC000120" / "msstats_input.tsv"
    ccle_p = REPO / "data" / "results" / "CCLE_corrected" / "msstats_input.tsv"
    if not cptac_p.exists() or not ccle_p.exists():
        print("Missing msstats_input", file=sys.stderr)
        return 1

    print("Reading CPTAC (capped)...")
    cn, cs, ln, ls = collect(cptac_p, args.max_rows)
    print("Reading CCLE (capped)...")
    cn2, cs2, ln2, ls2 = collect(ccle_p, args.max_rows)

    plot_four_panels(
        cn, cs, cn2, cs2, "Intensity", "Bridge vs sample — raw linear", "bridge_boxplots_raw", args.max_rows
    )
    plot_four_panels(
        ln, ls, ln2, ls2, "log2(intensity+1)", "Bridge vs sample — log2(+1)", "bridge_boxplots_log", args.max_rows
    )

    for ext in (".pdf", ".png"):
        src_log = OUT / f"bridge_boxplots_log{ext}"
        src_raw = OUT / f"bridge_boxplots_raw{ext}"
        if src_log.is_file() and src_raw.is_file():
            shutil.copy(src_log, OUT / f"sample_boxplots_log{ext}")
            shutil.copy(src_raw, OUT / f"sample_boxplots_raw{ext}")

    print("Wrote figures in", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
