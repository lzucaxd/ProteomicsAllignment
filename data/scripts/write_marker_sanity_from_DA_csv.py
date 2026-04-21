#!/usr/bin/env python3
"""Rebuild DA_MSstatsTMT_*_marker_sanity.csv from DA_MSstatsTMT_*.csv (same logic as R script)."""
import csv
import sys
from pathlib import Path

BASAL = {"KRT5", "KRT14", "KRT17", "EGFR", "FOXC1"}
LUM = {"ESR1", "GATA3", "FOXA1", "KRT18", "PGR"}
MARKERS = BASAL | LUM


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: write_marker_sanity_from_DA_csv.py <DA_MSstatsTMT_Luminal_vs_Basal.csv> <out_marker_sanity.csv>")
        sys.exit(1)
    inp = Path(sys.argv[1])
    outp = Path(sys.argv[2])
    rows_out = []
    with open(inp, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            g = (row.get("Gene_symbol") or "").strip()
            if g not in MARKERS:
                continue
            try:
                fc = float(row["log2FC"])
            except (KeyError, ValueError):
                continue
            try:
                ap = float(row.get("adj.pvalue") or row.get("adj.P.Val") or "")
            except ValueError:
                ap = float("nan")
            prot = row.get("Protein") or row.get("ProteinName") or ""
            exp = "Basal_up" if g in BASAL else "Luminal_up"
            # Contrast is Luminal - Basal: basal markers expect fc < 0, luminal fc > 0
            if fc > 0:
                obs = "Luminal_up"
            elif fc < 0:
                obs = "Basal_up"
            else:
                obs = "NS"
            rows_out.append(
                {
                    "Protein": prot,
                    "Gene_symbol": g,
                    "log2FC": fc,
                    "adj.pvalue": ap,
                    "expected": exp,
                    "observed": obs,
                    "direction_ok": str(exp == obs).upper(),
                }
            )
    rows_out.sort(key=lambda x: (x["Gene_symbol"], x["Protein"]))
    with open(outp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "Protein",
                "Gene_symbol",
                "log2FC",
                "adj.pvalue",
                "expected",
                "observed",
                "direction_ok",
            ],
        )
        w.writeheader()
        w.writerows(rows_out)
    print("Wrote", outp, "rows", len(rows_out))


if __name__ == "__main__":
    main()
