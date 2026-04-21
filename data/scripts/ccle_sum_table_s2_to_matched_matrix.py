#!/usr/bin/env python3
"""
Collapse CCLE Table S2 (Normalized Protein Expression) TenPx columns to CCLE codes
(median when multiple TenPx per code), then subset to MSI-mapped cell lines.

Input: glob data/ccle_sum/Table_S2*Protein*Normalized*.xlsx (first match)
Mapping: data/results/CCLE/msi_vs_mss/ccle_msi_label_mapping.csv (match_status == matched)

Output: data/results/CCLE/ccle_sum/table_s2_protein_matrix_cclecode_matched.csv.gz
        data/results/CCLE/ccle_sum/table_s2_matrix_build_summary.txt

Requires: pandas, openpyxl (e.g. pip install pandas openpyxl in project venv)

Run from repo root:
  .venv/bin/python data/scripts/ccle_sum_table_s2_to_matched_matrix.py
"""
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def find_xlsx() -> Path:
    for d in [ROOT / "data" / "ccle_sum", ROOT / "data" / "ccle_peptide"]:
        if not d.is_dir():
            continue
        for p in sorted(d.glob("Table_S2*Normalized*.xlsx")):
            return p
    raise FileNotFoundError("No Table_S2*Protein*Normalized*.xlsx under data/ccle_sum or data/ccle_peptide")


def main() -> None:
    xlsx = find_xlsx()
    map_path = ROOT / "data" / "results" / "CCLE" / "msi_vs_mss" / "ccle_msi_label_mapping.csv"
    out_dir = ROOT / "data" / "results" / "CCLE" / "ccle_sum"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_gz = out_dir / "table_s2_protein_matrix_cclecode_matched.csv.gz"

    with open(map_path, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r.get("match_status") == "matched" and r.get("msi_label") in ("MSI", "MSS")]
    wanted = {r["ccle_code_from_sample_info"].strip() for r in rows if r.get("ccle_code_from_sample_info")}

    meta = {
        "Protein_Id",
        "Gene_Symbol",
        "Description",
        "Group_ID",
        "Uniprot",
        "Uniprot_Acc",
    }
    print("Reading", xlsx, "...")
    df = pd.read_excel(xlsx, sheet_name="Normalized Protein Expression", engine="openpyxl")
    pep = [c for c in df.columns if str(c).endswith("_Peptides")]
    sample_cols = [c for c in df.columns if c not in meta and c not in pep]

    prefix = {c: re.sub(r"_TenPx\d+$", "", str(c)) for c in sample_cols}
    # collapse: median per prefix
    by_pref: dict[str, list[str]] = {}
    for c, p in prefix.items():
        by_pref.setdefault(p, []).append(c)

    collapsed = {}
    for p, cols in by_pref.items():
        sub = df[cols]
        collapsed[p] = sub.median(axis=1, skipna=True)

    mat = pd.DataFrame(collapsed)
    # subset to MSI workflow cell lines only
    have = [c for c in mat.columns if c in wanted]
    miss = sorted(wanted - set(have))
    mat = mat[have]

    meta_df = df[list(meta & set(df.columns))]
    out = pd.concat([meta_df, mat], axis=1)

    print("Writing", out_gz, "shape", out.shape)
    with gzip.open(out_gz, "wt", encoding="utf-8", newline="") as fh:
        out.to_csv(fh, index=False)

    summary = [
        "Table S2 → CCLE-code matrix (matched MSI/MSS labels only)",
        "=========================================================",
        f"Source xlsx: {xlsx.relative_to(ROOT)}",
        f"Rows (proteins): {len(out)}",
        f"Wanted CCLE codes (matched MSI/MSS): {len(wanted)}",
        f"Columns found in Table S2: {len(have)}",
        f"CCLE codes in mapping but missing from Table S2: {len(miss)}",
    ]
    if miss:
        summary.append("Missing codes (first 30):")
        summary.extend(f"  {m}" for m in miss[:30])
    summary.append("")
    summary.append("Collapse rule: columns named *CCLE*_TenPxNN → strip _TenPxNN; median if duplicates.")
    summary.append("Join key for DA: column name = ccle_code_from_sample_info in ccle_msi_label_mapping.csv")

    (out_dir / "table_s2_matrix_build_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
