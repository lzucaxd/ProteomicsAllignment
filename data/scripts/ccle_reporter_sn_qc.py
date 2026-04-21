#!/usr/bin/env python3
"""
QC summary of CCLE peptide-level reporter-ion S/N (rq_*_sn) from the input TSV.

Reads data/ccle_peptide/ccle_protein_quant_with_peptides_*.tsv in chunks.
Writes plots and tables under data/results/CCLE/qc_signal_to_noise/.

Usage (from repo root):
  python3 data/scripts/ccle_reporter_sn_qc.py

Or from data/:
  python3 scripts/ccle_reporter_sn_qc.py --tsv ccle_peptide/ccle_protein_quant_with_peptides_14745.tsv
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RQ_SN_PATTERN = re.compile(r"^rq_.*_sn$", re.IGNORECASE)
RUN_PLEX_PATTERN = re.compile(r"Prot_(\d+)", re.IGNORECASE)


def run_id_from_path(p: str) -> str:
    base = os.path.basename(p or "")
    return base.replace(".mzXML", "").replace(".raw", "").replace(".mzML", "").strip()


def mixture_from_run_id(run_id: str, valid_plex: set[int]) -> int | None:
    m = RUN_PLEX_PATTERN.search(run_id)
    if not m:
        return None
    mixture = int(m.group(1)) - 1
    if mixture in valid_plex:
        return mixture
    return None


def load_valid_plex_ids(sample_csv: Path) -> set[int]:
    out: set[int] = set()
    with sample_csv.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            key = "Protein 10-Plex ID"
            if key not in row:
                continue
            try:
                out.add(int((row[key] or "").strip()))
            except ValueError:
                continue
    return out


class Welford:
    __slots__ = ("n", "mean", "m2", "vmin", "vmax")

    def __init__(self) -> None:
        self.n = 0
        self.mean = 0.0
        self.m2 = 0.0
        self.vmin = float("inf")
        self.vmax = float("-inf")

    def update_many(self, x: np.ndarray) -> None:
        x = x[np.isfinite(x)]
        if x.size == 0:
            return
        for v in x.flat:
            self._update1(float(v))

    def _update1(self, v: float) -> None:
        self.n += 1
        d = v - self.mean
        self.mean += d / self.n
        d2 = v - self.mean
        self.m2 += d * d2
        if v < self.vmin:
            self.vmin = v
        if v > self.vmax:
            self.vmax = v

    def std(self) -> float:
        if self.n < 2:
            return float("nan")
        return math.sqrt(self.m2 / (self.n - 1))


def reservoir_update(
    reservoir: list[float],
    rng: random.Random,
    n_seen_finite: int,
    val: float,
    cap: int,
) -> None:
    """Algorithm R: uniform reservoir of size cap over a stream."""
    if len(reservoir) < cap:
        reservoir.append(val)
        return
    j = rng.randint(0, n_seen_finite - 1)
    if j < cap:
        reservoir[j] = val


def main() -> int:
    ap = argparse.ArgumentParser(description="CCLE rq_*_sn QC")
    root = Path(__file__).resolve().parents[2]
    default_tsv = root / "data" / "ccle_peptide" / "ccle_protein_quant_with_peptides_14745.tsv"
    default_sample = root / "data" / "ccle_peptide" / "sample_info_ccle.csv"
    default_out = root / "data" / "results" / "CCLE" / "qc_signal_to_noise"
    ap.add_argument("--tsv", type=Path, default=default_tsv)
    ap.add_argument("--sample_csv", type=Path, default=default_sample)
    ap.add_argument("--outdir", type=Path, default=default_out)
    ap.add_argument("--chunksize", type=int, default=200_000)
    ap.add_argument("--reservoir", type=int, default=250_000)
    ap.add_argument("--hist-bins", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if not args.tsv.is_file():
        print(f"Missing input TSV: {args.tsv}", file=sys.stderr)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    valid_plex = load_valid_plex_ids(args.sample_csv) if args.sample_csv.is_file() else set()

    head = pd.read_csv(args.tsv, sep="\t", nrows=0)
    sn_cols = [c for c in head.columns if RQ_SN_PATTERN.match(str(c).strip())]
    if not sn_cols:
        print("No rq_*_sn columns found.", file=sys.stderr)
        return 1

    usecols = sn_cols + ["RunLoadPath"]

    global_w = Welford()
    by_channel = {c: Welford() for c in sn_cols}
    by_run: dict[str, Welford] = defaultdict(Welford)
    by_mix: dict[int, Welford] = defaultdict(Welford)

    reservoir: list[float] = []
    n_seen_finite = 0
    total_cells = 0
    n_nan = 0

    chunk_iter = pd.read_csv(
        args.tsv,
        sep="\t",
        usecols=lambda c: c in set(usecols),
        chunksize=args.chunksize,
        low_memory=False,
    )

    for chunk in chunk_iter:
        chunk["run"] = chunk["RunLoadPath"].astype(str).map(run_id_from_path)

        for col in sn_cols:
            v = pd.to_numeric(chunk[col], errors="coerce").to_numpy()
            total_cells += v.size
            n_nan += int(np.isnan(v).sum())
            global_w.update_many(v)
            by_channel[col].update_many(v)
            for val in v[np.isfinite(v)].flat:
                n_seen_finite += 1
                reservoir_update(reservoir, rng, n_seen_finite, float(val), args.reservoir)

        # Per-run and mixture (vectorized groupby)
        for rid, sub in chunk.groupby("run", sort=False):
            stacked = pd.to_numeric(sub[sn_cols].stack(), errors="coerce").values
            stacked = stacked[np.isfinite(stacked)]
            if stacked.size == 0:
                continue
            by_run[rid].update_many(stacked)
            mix = mixture_from_run_id(rid, valid_plex) if valid_plex else None
            if mix is not None:
                by_mix[mix].update_many(stacked)

    # Histogram from reservoir (defensible global view; full table is too large to load at once)
    rs = np.array(reservoir, dtype=float)
    rs.sort()
    p_lo, p_hi = np.percentile(rs, [0.5, 99.5])
    edges = np.linspace(p_lo, p_hi, args.hist_bins + 1)
    hist_counts, _ = np.histogram(rs, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(centers, hist_counts, width=widths, align="center", alpha=0.8, color="steelblue", edgecolor="none")
    ax.set_xlabel("Reporter-ion S/N (rq_*_sn)")
    ax.set_ylabel("Count")
    ax.set_title(
        f"CCLE reporter-ion S/N (uniform reservoir n={len(reservoir)}; axis trimmed p0.5–p99.5)"
    )
    ax2 = ax.twinx()
    dens, _ = np.histogram(rs, bins=edges, density=True)
    ax2.plot(centers, dens, color="darkorange", lw=1.8, alpha=0.95, label="Density")
    ax2.set_ylabel("Density")
    ax2.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.outdir / "qc_sn_histogram_density.png", dpi=150)
    plt.close(fig)

    # Summary tables
    rows = [
        {
            "scope": "overall_all_channels_pooled",
            "n_values": global_w.n,
            "mean": global_w.mean,
            "std": global_w.std(),
            "min": global_w.vmin if global_w.n else float("nan"),
            "p01": float(np.percentile(rs, 1)) if rs.size > 50 else float("nan"),
            "p50": float(np.percentile(rs, 50)) if rs.size > 50 else float("nan"),
            "p99": float(np.percentile(rs, 99)) if rs.size > 50 else float("nan"),
            "max": global_w.vmax if global_w.n else float("nan"),
            "frac_nan": n_nan / total_cells if total_cells else float("nan"),
        }
    ]
    for col in sn_cols:
        w = by_channel[col]
        rows.append(
            {
                "scope": f"channel:{col}",
                "n_values": w.n,
                "mean": w.mean,
                "std": w.std(),
                "min": w.vmin if w.n else float("nan"),
                "p01": float("nan"),
                "p50": float("nan"),
                "p99": float("nan"),
                "max": w.vmax if w.n else float("nan"),
                "frac_nan": float("nan"),
            }
        )
    pd.DataFrame(rows).to_csv(args.outdir / "qc_sn_summary.tsv", sep="\t", index=False)

    ch_rows = []
    for col in sn_cols:
        w = by_channel[col]
        ch_rows.append(
            {
                "channel_col": col,
                "n_values": w.n,
                "mean": w.mean,
                "std": w.std(),
                "min": w.vmin,
                "max": w.vmax,
            }
        )
    pd.DataFrame(ch_rows).to_csv(args.outdir / "qc_sn_by_channel.tsv", sep="\t", index=False)

    run_rows = []
    run_means: list[float] = []
    for rid, w in sorted(by_run.items(), key=lambda x: x[0]):
        run_means.append(w.mean)
        run_rows.append(
            {
                "run_id": rid,
                "mixture": mixture_from_run_id(rid, valid_plex) if valid_plex else "",
                "n_values": w.n,
                "mean_sn": w.mean,
                "std_sn": w.std(),
                "min_sn": w.vmin,
                "max_sn": w.vmax,
            }
        )
    pd.DataFrame(run_rows).to_csv(args.outdir / "qc_sn_by_run.tsv", sep="\t", index=False)

    mx_rows = []
    for mix in sorted(by_mix.keys()):
        w = by_mix[mix]
        mx_rows.append(
            {
                "mixture": mix,
                "n_values": w.n,
                "mean_sn": w.mean,
                "std_sn": w.std(),
                "min_sn": w.vmin,
                "max_sn": w.vmax,
            }
        )
    pd.DataFrame(mx_rows).to_csv(args.outdir / "qc_sn_by_mixture.tsv", sep="\t", index=False)

    rm = np.array(run_means, dtype=float)
    if rm.size:
        fig2, axb = plt.subplots(figsize=(8, 3))
        axb.boxplot(rm, vert=True)
        axb.set_ylabel("Mean S/N per run (pooled channels)")
        axb.set_title("Distribution of run-level mean reporter-ion S/N")
        fig2.tight_layout()
        fig2.savefig(args.outdir / "qc_sn_run_means_boxplot.png", dpi=150)
        plt.close(fig2)

    overall_mean = global_w.mean
    note = f"""CCLE reporter-ion S/N QC (rq_*_sn columns)
============================================
Input: {args.tsv.name}
Approx. table rows (cells / n_channels): {total_cells // max(len(sn_cols), 1)}.

What was computed
-----------------
- Pooled all finite rq_*_sn values across channels and rows (one value = one reporter measurement).
- Global mean/std/min/max over all finite values (n={global_w.n}).
- Approximate p01/p50/p99 from a uniform random reservoir (n={len(reservoir)}) for speed on large files.
- Histogram/density plot uses the same reservoir; x-axis trimmed to p0.5–p99.5 of the reservoir so the figure
  shows the bulk of the distribution (extremes still in min/max and quantile table).
- Per-channel: mean/std/min/max over all rows (no reservoir).
- Per-run: run_id = basename(RunLoadPath); stats pool all channels and rows for that run.
- Mixture: int(NN)-1 from 'Prot_NN' in the run name when it matches Protein 10-Plex IDs in sample_info_ccle.csv.

Model / units
-------------
Values are instrument-reported reporter ion signal-to-noise as provided in the CCLE export (not re-derived here).

Reasonableness
--------------
- Overall mean S/N ≈ {overall_mean:.4g} (finite values n={global_w.n}).
- Missing/non-finite fraction ≈ {n_nan} / {total_cells} ({100 * n_nan / max(total_cells, 1):.3f}%).
- Runs summarized: {len(by_run)}; mixtures with mapped Prot_NN: {len(by_mix)}.

Red flags vs reassuring signs
-----------------------------
- A left tail of low S/N is expected (weak peptides); compare p01 to median.
- If several runs have mean S/N far below the cohort (qc_sn_by_run.tsv), inspect those LC-MS files.
- Channel means should be broadly comparable; large systematic shifts may indicate batch/plex effects upstream.

Outputs
-------
- qc_sn_summary.tsv, qc_sn_by_channel.tsv, qc_sn_by_run.tsv, qc_sn_by_mixture.tsv
- qc_sn_histogram_density.png, qc_sn_run_means_boxplot.png
"""
    (args.outdir / "qc_sn_interpretation.txt").write_text(note, encoding="utf-8")

    print(f"Wrote QC under {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
