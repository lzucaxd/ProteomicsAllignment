#!/usr/bin/env python3
"""
Build effect_size_and_spread_summary.tsv: gene-level log2FC, adj.P, SE (uncertainty proxy)
from collapsed protein-level DA rows (same rule as build_benchmark_v1_artifacts.py).

SE column in MSstatsTMT output = standard error of log2FC (per protein).

Usage:
  python3 data/scripts/build_effect_size_spread_summary.py
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "reports" / "benchmark_v1" / "diagnostics_feedback"

CPTAC_DA = REPO / "data" / "results" / "PDC000120" / "DA_subtype_subset_runs" / "DA_MSstatsTMT_Luminal_vs_Basal.csv"
CCLE_DA = REPO / "data" / "results" / "CCLE_corrected" / "DA_luminal_vs_basal" / "DA_MSstatsTMT_Luminal_vs_Basal.csv"
SHARED = REPO / "reports" / "benchmark_v1" / "shared_feature_table.csv"

MARKERS = {
    "ESR1",
    "GATA3",
    "FOXA1",
    "PGR",
    "KRT18",
    "KRT5",
    "KRT14",
    "KRT17",
    "EGFR",
    "FOXC1",
}


def load_best(path: Path) -> dict:
    rows_by_gene: dict[str, list[dict]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f):
            if (row.get("Protein") or "").startswith("##"):
                continue
            g = (row.get("Gene_symbol") or "").strip()
            if not g:
                continue
            try:
                rows_by_gene[g].append(
                    {
                        "log2FC": float(row["log2FC"]),
                        "adj": float(row["adj.pvalue"]),
                        "SE": float(row.get("SE") or "nan"),
                    }
                )
            except (ValueError, KeyError):
                continue
    best = {}
    for g, lst in rows_by_gene.items():

        def k(d):
            return (d["adj"], d["SE"] if d["SE"] == d["SE"] else 1e9, -abs(d["log2FC"]))

        lst.sort(key=k)
        b = lst[0]
        best[g] = b
    return best


def main() -> int:
    cptac = load_best(CPTAC_DA)
    ccle = load_best(CCLE_DA)
    genes = []
    with SHARED.open(encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            genes.append(row["GeneSymbol"])

    rows = []
    for g in genes:
        ca, cb = cptac.get(g), ccle.get(g)
        if not ca or not cb:
            continue
        same = (ca["log2FC"] > 0 and cb["log2FC"] > 0) or (ca["log2FC"] < 0 and cb["log2FC"] < 0)
        if ca["log2FC"] == 0 or cb["log2FC"] == 0:
            same = False
        rows.append(
            {
                "GeneSymbol": g,
                "CPTAC_log2FC": f"{ca['log2FC']:.10g}",
                "CCLE_log2FC": f"{cb['log2FC']:.10g}",
                "CPTAC_adj_pvalue": f"{ca['adj']:.10g}",
                "CCLE_adj_pvalue": f"{cb['adj']:.10g}",
                "CPTAC_significant_FDR05": "yes" if ca["adj"] < 0.05 else "no",
                "CCLE_significant_FDR05": "yes" if cb["adj"] < 0.05 else "no",
                "same_direction": "yes" if same else "no",
                "CPTAC_SE_log2FC": f"{ca['SE']:.10g}" if ca["SE"] == ca["SE"] else "",
                "CCLE_SE_log2FC": f"{cb['SE']:.10g}" if cb["SE"] == cb["SE"] else "",
                "marker_flag": "yes" if g in MARKERS else "no",
                "spread_proxy_note": "SE is MSstatsTMT-reported SE of log2FC per collapsed protein row",
            }
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "effect_size_and_spread_summary.tsv"
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        w.writerows(rows)
    print("Wrote", out, "rows", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
