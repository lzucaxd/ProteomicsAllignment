#!/usr/bin/env python3
"""
Export the second sheet of Table_S1_Sample_Information (CCLE) to CSV.
Uses only stdlib (zipfile + xml). Run from data/ccle_peptide/:
  python3 export_sample_sheet2_csv.py
  -> writes sample_info_ccle.csv
"""
import csv
import os
import re
import xml.etree.ElementTree as ET
import zipfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(SCRIPT_DIR, "Table_S1_Sample_Information (1).xlsx")
OUT_CSV = os.path.join(SCRIPT_DIR, "sample_info_ccle.csv")


def read_shared_strings(zipf):
    with zipf.open("xl/sharedStrings.xml") as f:
        root = ET.parse(f).getroot()
    uri = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    strings = []
    for si in root.findall(f"{{{uri}}}si"):
        t = si.find(f"{{{uri}}}t")
        if t is not None and t.text:
            strings.append(t.text)
        else:
            r = si.find(f"{{{uri}}}r")
            if r is not None:
                parts = []
                for t in r.findall(f"{{{uri}}}t"):
                    if t.text:
                        parts.append(t.text)
                strings.append("".join(parts))
            else:
                strings.append("")
    return strings


def cell_ref_to_col(cr):
    """A1 -> 0, B1 -> 1, AA1 -> 26 (0-based)."""
    m = re.match(r"^([A-Z]+)", cr.upper())
    if not m:
        return 0
    col = 0
    for c in m.group(1):
        col = col * 26 + (ord(c) - ord("A"))
    return col


def read_sheet(zipf, sheet_path, shared_strings):
    with zipf.open(sheet_path) as f:
        root = ET.parse(f).getroot()
    # Excel uses default namespace; tags are {uri}localname
    uri = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    sheet_data = root.find(f".//{{{uri}}}sheetData")
    if sheet_data is None:
        return []
    rows = []
    for row in sheet_data.findall(f"{{{uri}}}row"):
        r = row.get("r")
        cells = {}
        for c in row.findall(f"{{{uri}}}c"):
            ref = c.get("r", "")
            col = cell_ref_to_col(ref)
            t = c.get("t")
            v = c.find(f"{{{uri}}}v")
            val = v.text if v is not None and v.text else ""
            if t == "s" and val.isdigit():
                idx = int(val)
                val = shared_strings[idx] if idx < len(shared_strings) else ""
            cells[col] = val
        max_col = max(cells.keys()) if cells else -1
        row_list = [cells.get(j, "") for j in range(max_col + 1)]
        rows.append(row_list)
    return rows


def main():
    if not os.path.isfile(XLSX):
        raise FileNotFoundError(XLSX)
    with zipfile.ZipFile(XLSX, "r") as z:
        strings = read_shared_strings(z)
        rows = read_sheet(z, "xl/worksheets/sheet2.xml", strings)
    if not rows:
        raise SystemExit("Sheet2 is empty.")
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        for row in rows:
            w.writerow(["" if c is None else str(c).strip() for c in row])
    print("Wrote", OUT_CSV, "(", len(rows), "rows)")


if __name__ == "__main__":
    main()
