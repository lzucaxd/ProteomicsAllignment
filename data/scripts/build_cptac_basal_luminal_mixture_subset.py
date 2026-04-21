#!/usr/bin/env python3
"""
Build sample-level mixture balance for CPTAC PDC000120 Basal vs Luminal (LumA+LumB),
subset annotations, and companion CSVs for mixture-balanced DA.

Rule for keep_for_subset (documented in outputs):
  - DROP if either subtype count is 0 (one subtype absent from the TMT mixture).
  - DROP if min(n_Basal, n_Luminal) == 1 AND total >= 6 (singleton from one subtype in
    a 6+ sample mixture — "nearly absent" minority arm for stable within-plex comparison).

This does not change the biological contrast (Basal vs Luminal); it restricts which
bioreplicates enter the analysis so mixture is not perfectly aligned with subtype
in problematic plexes.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


def luminal_group(pam50: str) -> str | None:
    p = (pam50 or "").strip()
    if p in ("LumA", "LumB"):
        return "Luminal"
    if p == "Basal":
        return "Basal"
    return None


def decide_keep(n_b: int, n_l: int) -> tuple[str, str]:
    total = n_b + n_l
    if n_b < 1 or n_l < 1:
        return "drop", "missing_subtype"
    if min(n_b, n_l) == 1 and total >= 6:
        return "drop", "near_absent_minority"
    return "keep", ""


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    res = repo / "data" / "results" / "PDC000120"
    tumor_path = res / "DA_subtype_tumor_only.csv"
    sample_path = res / "DA_sample_annotation.csv"
    out_balance = res / "cptac_basal_luminal_mixture_balance.csv"
    out_keepdrop = res / "cptac_basal_luminal_subset_keepdrop.csv"
    out_tumor_subset = res / "DA_subtype_tumor_only_basal_luminal_subset.csv"
    out_sample_subset = res / "DA_sample_annotation_basal_luminal_subset.csv"

    rows = list(csv.DictReader(open(tumor_path, newline="", encoding="utf-8")))
    mix_counts: dict[str, dict[str, int]] = defaultdict(lambda: {"Basal": 0, "Luminal": 0})
    for r in rows:
        g = luminal_group(r.get("pam50") or "")
        if g is None:
            continue
        m = (r.get("mixture") or "").strip()
        if not m:
            continue
        mix_counts[m][g] += 1

    mixtures_sorted = sorted(mix_counts.keys())
    balance_rows = []
    keepdrop_rows = []
    for m in mixtures_sorted:
        c = mix_counts[m]
        nb, nl = c["Basal"], c["Luminal"]
        tot = nb + nl
        balance_rows.append({"mixture": m, "n_Basal": nb, "n_Luminal": nl, "total": tot})
        k, reason = decide_keep(nb, nl)
        keepdrop_rows.append(
            {
                "mixture": m,
                "n_Basal": nb,
                "n_Luminal": nl,
                "total": tot,
                "keep_for_subset": k,
                "drop_reason": reason or "none",
            }
        )

    with open(out_balance, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["mixture", "n_Basal", "n_Luminal", "total"])
        w.writeheader()
        w.writerows(balance_rows)

    with open(out_keepdrop, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "mixture",
                "n_Basal",
                "n_Luminal",
                "total",
                "keep_for_subset",
                "drop_reason",
            ],
        )
        w.writeheader()
        w.writerows(keepdrop_rows)

    kept_mix = {r["mixture"] for r in keepdrop_rows if r["keep_for_subset"] == "keep"}

    def filter_table(path: Path, outp: Path) -> int:
        rdr = csv.DictReader(open(path, newline="", encoding="utf-8"))
        fieldnames = rdr.fieldnames
        if not fieldnames:
            raise SystemExit(f"No header: {path}")
        out_lines = []
        for r in rdr:
            pam = (r.get("pam50") or "").strip()
            g = luminal_group(pam)
            if g is None:
                continue
            m = (r.get("mixture") or "").strip()
            if m in kept_mix:
                out_lines.append(r)
        with open(outp, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(out_lines)
        return len(out_lines)

    n_tumor = filter_table(tumor_path, out_tumor_subset)
    n_sample = filter_table(sample_path, out_sample_subset)

    print("Wrote:", out_balance)
    print("Wrote:", out_keepdrop)
    print("Wrote:", out_tumor_subset, f"({n_tumor} rows)")
    print("Wrote:", out_sample_subset, f"({n_sample} rows)")
    print("Kept mixtures:", len(kept_mix), "/", len(mixtures_sorted))


if __name__ == "__main__":
    main()
