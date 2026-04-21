#!/usr/bin/env python3
"""
Build CCLE lineage contrast candidate table and conservative choice file.
Uses data/ccle_peptide/sample_info_ccle.csv (Tissue of Origin) and
data/results/CCLE/gene_matrix.csv column names.

Run from repo root:
  python3 data/scripts/ccle_lineage_contrast_plan.py
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GM = ROOT / "data/results/CCLE/gene_matrix.csv"
SAMPLE_INFO = ROOT / "data/ccle_peptide/sample_info_ccle.csv"
OUT_CSV = ROOT / "data/results/CCLE/ccle_lineage_contrast_candidates.csv"
OUT_COUNTS = ROOT / "data/results/CCLE/ccle_lineage_counts_by_bucket.csv"
OUT_TXT = ROOT / "data/results/CCLE/ccle_lineage_contrast_choice.txt"


def lineage_bucket(tissue: str) -> str | None:
    t = tissue.strip()
    m = {
        "Breast": "breast",
        "Kidney": "kidney",
        "Lung": "lung",
        "Pancreas": "pancreas",
        "Stomach": "stomach",
        "Central Nervous System": "cns",
    }
    return m.get(t)


def main() -> None:
    with open(GM) as f:
        r = csv.reader(f)
        header = next(r)
        cols = header[2:]
        na = {c: 0 for c in cols}
        n_genes = 0
        for row in r:
            n_genes += 1
            for i, c in enumerate(cols):
                v = row[i + 2].strip()
                if v == "":
                    na[c] += 1
                else:
                    try:
                        float(v)
                    except ValueError:
                        na[c] += 1
    na_frac = {c: na[c] / max(n_genes, 1) for c in cols}

    # cell line -> tissue, mixture (from sample sheet; first row wins if duplicate)
    info: dict[str, dict] = {}
    with open(SAMPLE_INFO) as f:
        r = csv.DictReader(f)
        for row in r:
            cl = row["Cell Line"]
            if cl not in info:
                info[cl] = {
                    "tissue": row["Tissue of Origin"],
                    "mixture": row["Protein 10-Plex ID"].strip(),
                }

    line_of = {}
    for c in cols:
        if c not in info:
            raise SystemExit(f"Column {c!r} not in sample_info")
        b = lineage_bucket(info[c]["tissue"])
        if b is None:
            continue
        line_of[c] = b

    def line_cols(bucket: str) -> list[str]:
        return [c for c in cols if line_of.get(c) == bucket]

    def mixture_set(bucket: str) -> set[str]:
        return {info[c]["mixture"] for c in line_cols(bucket)}

    def mean_na(bucket: str) -> float:
        lc = line_cols(bucket)
        if not lc:
            return float("nan")
        return sum(na_frac[c] for c in lc) / len(lc)

    contrasts = [
        ("breast_vs_kidney", "breast", "kidney"),
        ("breast_vs_lung", "breast", "lung"),
        ("pancreas_vs_stomach", "pancreas", "stomach"),
        ("cns_vs_lung", "cns", "lung"),
    ]

    rows_out = []
    for cid, a, b in contrasts:
        ca, cb = line_cols(a), line_cols(b)
        sa, sb = mixture_set(a), mixture_set(b)
        inter = len(sa & sb)
        union = len(sa | sb)
        jacc = inter / union if union else 0.0
        n_a, n_b = len(ca), len(cb)
        ratio = max(n_a, n_b) / max(min(n_a, n_b), 1)

        # Heuristic score 1–5
        score = 3
        notes = []
        if n_a < 5 or n_b < 5:
            score = min(score, 1)
            notes.append("very_small_group")
        if ratio > 4:
            score = min(score, 2)
            notes.append("strong_n_imbalance")
        elif ratio > 2.5:
            score = min(score, 3)
            notes.append("moderate_n_imbalance")
        if jacc < 0.2:
            score = min(score, 2)
            notes.append("low_shared_TMT_mixtures_batch_risk")
        elif jacc >= 0.45:
            notes.append("good_mixture_overlap_across_groups")
            score = max(score, 4)
        if jacc >= 0.2 and jacc < 0.45:
            notes.append("partial_mixture_overlap")

        rec = (
            "Reasonable first-pass benchmark if interpreted carefully."
            if score >= 3
            else "Weak for batch-safe lineage DA; use heavy caveats or avoid."
        )
        if cid == "breast_vs_lung" and score >= 4:
            rec = "strongest among four on TMT mixture overlap + adequate n; n imbalance remains."
        if cid == "breast_vs_kidney":
            rec = (
                "Biologically clean contrast but low shared mixtures with lung/breast cohorts "
                "implies batch structure may dominate unless modeled."
            )
        if cid == "cns_vs_lung":
            rec = "Very imbalanced n (13 vs 76); not recommended as primary."

        rows_out.append(
            {
                "contrast_id": cid,
                "group_a": a,
                "group_b": b,
                "n_samples_a": n_a,
                "n_samples_b": n_b,
                "n_unique_cell_lines_a": n_a,
                "n_unique_cell_lines_b": n_b,
                "mean_na_frac_a": round(mean_na(a), 5),
                "mean_na_frac_b": round(mean_na(b), 5),
                "n_mixtures_a": len(sa),
                "n_mixtures_b": len(sb),
                "n_mixtures_shared": inter,
                "mixture_jaccard": round(jacc, 5),
                "n_imbalance_ratio_max_min": round(ratio, 3),
                "imbalance_note": ";".join(notes) if notes else "none",
                "recommendation_score_1to5": score,
                "recommendation_note": rec,
            }
        )

    # Per-lineage counts (matrix columns only)
    buckets = ["breast", "kidney", "lung", "pancreas", "stomach", "cns"]
    count_rows = []
    for b in buckets:
        lc = line_cols(b)
        count_rows.append(
            {
                "lineage_bucket": b,
                "n_cell_lines_in_gene_matrix": len(lc),
                "notes": "Tissue from sample_info_ccle.csv Tissue of Origin",
            }
        )
    with open(OUT_COUNTS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(count_rows[0].keys()))
        w.writeheader()
        w.writerows(count_rows)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)

    chosen = "breast_vs_lung"
    rationale = """
Chosen contrast (conservative): breast_vs_lung (Breast vs Lung)

Why not breast_vs_kidney (user preference when feasible):
  - TMT mixture overlap between breast and kidney is very low (see CSV; Jaccard typically ~0.17).
  - That pattern means breast and kidney lines largely occupy different MS/TMT experiment
    blocks; lineage differences would be strongly entangled with batch / run structure
    in a simple limma contrast without joint modeling.

Why breast_vs_lung:
  - Largest overlap of TMT mixtures shared between the two lineages (Jaccard typically ~0.50),
    so both tissues are represented across many of the same plexes — better for a
    batch-aware exploratory benchmark than kidney vs breast here.
  - Sample counts: ~30 breast vs ~76 lung lines (imbalanced but both large; imbalance
    is explicit in interpretation).

Caveats:
  - This is still cell-line proteomics with one column per line; residual batch/TMT
    effects are possible.
  - Do not over-interpret significance; use effect sizes and known tissue markers.

Alternatives:
  - pancreas_vs_stomach: more balanced n (~19 vs ~14) but lower mixture overlap than
    breast_vs_lung; reasonable secondary benchmark.
  - cns_vs_lung: not recommended due to extreme n imbalance.

Metadata source: data/ccle_peptide/sample_info_ccle.csv column "Tissue of Origin".
"""

    OUT_TXT.write_text(rationale.strip() + "\n")
    print("Wrote", OUT_COUNTS)
    print("Wrote", OUT_CSV)
    print("Wrote", OUT_TXT)


if __name__ == "__main__":
    main()
