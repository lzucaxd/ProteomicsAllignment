#!/usr/bin/env python3
"""
Build CCLE proteomics ↔ DepMap MSI labels (MSI if MSIScore > 20, else MSS).

Inputs (auto-detected under repo root unless overridden):
  - DepMap Model.csv (expects columns including ModelID, CCLEName, OncotreeLineage, ...)
  - DepMap OmicsGlobalSignatures.csv (MSIScore; uses rows with IsDefaultEntryForModel == Yes)
  - CCLE sample_info_ccle.csv (Cell Line, CCLE Code)
  - CCLE gene_matrix.csv (column names = cell line names)

Matching rules (conservative):
  1. Matrix column name must equal "Cell Line" in sample_info (exact string match).
  2. sample_info "CCLE Code" must equal DepMap Model "CCLEName" (exact match).
  3. Join DepMap signatures on ModelID; require IsDefaultEntryForModel == Yes for MSIScore.
  4. If Model exists but ModelID has no default signature row → MSIScore missing (no label).
  5. If no Model row for CCLEName → unmatched (should not occur if CCLE internal consistency holds).

Run from repo root:
  python3 data/scripts/build_ccle_msi_labels.py

Optional:
  python3 data/scripts/build_ccle_msi_labels.py --repo-root /path/to/ProteomicsAllignment
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


def find_depmap_files(root: Path) -> tuple[Path, Path]:
    """Locate Model.csv and OmicsGlobalSignatures.csv under common folder names."""
    candidates = [
        root / "data" / "MSI VS MSS",
        root / "data" / "msi_vs_mss",
        root / "msi_vs_mss",
    ]
    model = sig = None
    for base in candidates:
        if not base.is_dir():
            continue
        for p in base.iterdir():
            if p.is_file() and re.match(r"(?i)model\.csv$", p.name):
                model = p
            if p.is_file() and "omics" in p.name.lower() and "signature" in p.name.lower() and p.suffix.lower() == ".csv":
                sig = p
        if model and sig:
            return model, sig
    raise FileNotFoundError(
        "Could not find DepMap Model.csv and OmicsGlobalSignatures.csv. "
        "Place them under data/MSI VS MSS/ or data/msi_vs_mss/."
    )


def load_matrix_columns(gm_path: Path) -> list[str]:
    with open(gm_path, newline="") as f:
        row = next(csv.reader(f))
    return row[2:]


def main() -> None:
    ap = argparse.ArgumentParser(description="Build CCLE MSI label table from DepMap + CCLE metadata.")
    ap.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Repository root (default: parent of data/)",
    )
    ap.add_argument(
        "--msi-threshold",
        type=float,
        default=20.0,
        help="MSI if MSIScore > this value (default 20)",
    )
    args = ap.parse_args()
    root = args.repo_root.resolve()

    model_path, sig_path = find_depmap_files(root)
    sample_info = root / "data" / "ccle_peptide" / "sample_info_ccle.csv"
    gene_matrix = root / "data" / "results" / "CCLE" / "gene_matrix.csv"
    out_dir = root / "data" / "results" / "CCLE" / "msi_vs_mss"
    out_dir.mkdir(parents=True, exist_ok=True)

    schema_path = out_dir / "depmap_input_schema.txt"
    with open(model_path, newline="") as f:
        model_rows = list(csv.DictReader(f))
    with open(sig_path, newline="") as f:
        sig_reader = csv.DictReader(f)
        sig_fieldnames = sig_reader.fieldnames or []
        sig_rows = list(sig_reader)

    # Schema notes (compact)
    model_keys = list(model_rows[0].keys()) if model_rows else []
    id_like = [k for k in model_keys if k in ("ModelID", "CCLEName", "StrippedCellLineName", "CellLineName", "PatientID")]
    lineage_like = [k for k in model_keys if "Oncotree" in k or "Lineage" in k]
    sex_like = [k for k in model_keys if k.lower() == "sex"]
    schema_lines = [
        "DepMap + CCLE — compact schema summary",
        "=====================================",
        f"Model file: {model_path.relative_to(root)}",
        f"  All columns ({len(model_keys)}): {', '.join(model_keys)}",
        f"  Likely ID columns: {', '.join(id_like) if id_like else '(none detected)'}",
        f"  Likely lineage / disease columns: {', '.join(lineage_like)}",
        f"  Sex column: {sex_like[0] if sex_like else 'not found'}",
        "",
        f"Signatures file: {sig_path.relative_to(root)}",
        f"  All columns: {', '.join(sig_fieldnames)}",
        "  Key: ModelID (join to Model.ModelID), MSIScore (omics MSI signature score),",
        "       IsDefaultEntryForModel (use Yes only for one row per model).",
        "  Note: ProfileID is not used here; SequencingID is per omics profile.",
        "",
        "CCLE inputs",
        "===========",
        f"sample_info: {sample_info.relative_to(root)} — Cell Line (display name), CCLE Code (DepMap CCLEName key)",
        f"gene_matrix: {gene_matrix.relative_to(root)} — column names = Cell Line",
        "",
        "Intermediate mapping: none required; join path is",
        "  gene_matrix column → sample_info.Cell Line → sample_info.CCLE Code → Model.CCLEName.",
        "",
    ]
    schema_path.write_text("\n".join(schema_lines) + "\n")

    # Model by CCLEName (unique)
    by_ccle: dict[str, dict] = {}
    for r in model_rows:
        cn = (r.get("CCLEName") or "").strip()
        if not cn:
            continue
        if cn in by_ccle:
            raise SystemExit(f"Duplicate CCLEName in Model.csv: {cn}")
        by_ccle[cn] = r

    # Signature: default entry per ModelID
    sig_by_model: dict[str, dict] = {}
    for r in sig_rows:
        if r.get("IsDefaultEntryForModel") != "Yes":
            continue
        mid = r.get("ModelID", "").strip()
        if not mid:
            continue
        if mid in sig_by_model:
            raise SystemExit(f"Multiple default signature rows for ModelID {mid}")
        sig_by_model[mid] = r

    # sample_info: first occurrence per Cell Line (same as CCLE pipelines)
    cl_to_code: dict[str, str] = {}
    with open(sample_info, newline="") as f:
        for r in csv.DictReader(f):
            cl = (r.get("Cell Line") or "").strip()
            if cl and cl not in cl_to_code:
                cl_to_code[cl] = (r.get("CCLE Code") or "").strip()

    matrix_cols = load_matrix_columns(gene_matrix)
    msi_threshold = args.msi_threshold

    mapping_rows: list[dict] = []
    for col in matrix_cols:
        row: dict = {
            "proteomics_column_name": col,
            "ccle_code_from_sample_info": "",
            "match_status": "",
            "depmap_model_id": "",
            "depmap_cell_line_name": "",
            "depmap_stripped_cell_line_name": "",
            "depmap_oncotree_lineage": "",
            "depmap_oncotree_primary_disease": "",
            "depmap_sex": "",
            "MSIScore": "",
            "msi_label": "",
            "notes": "",
        }
        code = cl_to_code.get(col)
        if not code:
            row["match_status"] = "no_sample_info_for_column"
            row["notes"] = "Cell line name not found in sample_info_ccle.csv"
            mapping_rows.append(row)
            continue
        row["ccle_code_from_sample_info"] = code
        m = by_ccle.get(code)
        if not m:
            row["match_status"] = "sample_info_ok_no_depmap_model"
            row["notes"] = f"No Model row with CCLEName == {code!r}"
            mapping_rows.append(row)
            continue
        mid = m.get("ModelID", "").strip()
        row["depmap_model_id"] = mid
        row["depmap_cell_line_name"] = m.get("CellLineName", "")
        row["depmap_stripped_cell_line_name"] = m.get("StrippedCellLineName", "")
        row["depmap_oncotree_lineage"] = m.get("OncotreeLineage", "")
        row["depmap_oncotree_primary_disease"] = m.get("OncotreePrimaryDisease", "")
        row["depmap_sex"] = (m.get("Sex") or "").strip()
        sr = sig_by_model.get(mid)
        if not sr:
            row["match_status"] = "depmap_model_ok_no_default_msi_signature"
            row["notes"] = "ModelID not in OmicsGlobalSignatures with IsDefaultEntryForModel==Yes"
            mapping_rows.append(row)
            continue
        raw = (sr.get("MSIScore") or "").strip()
        if raw == "":
            row["match_status"] = "signature_row_missing_MSIScore"
            mapping_rows.append(row)
            continue
        score = float(raw)
        row["MSIScore"] = f"{score:.6g}"
        row["msi_label"] = "MSI" if score > msi_threshold else "MSS"
        row["match_status"] = "matched"
        mapping_rows.append(row)

    map_path = out_dir / "ccle_msi_label_mapping.csv"
    with open(map_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mapping_rows[0].keys()))
        w.writeheader()
        w.writerows(mapping_rows)

    matched = [r for r in mapping_rows if r["match_status"] == "matched"]
    labels = [r["msi_label"] for r in matched]
    ctr = Counter(labels)
    lineage_ctr = Counter(r["depmap_oncotree_lineage"] for r in matched)
    lineage_msi = Counter()
    for r in matched:
        lineage_msi[f'{r["depmap_oncotree_lineage"]}|{r["msi_label"]}'] += 1

    counts_path = out_dir / "ccle_msi_group_counts.txt"
    counts_lines = [
        "CCLE proteomics columns: %d" % len(matrix_cols),
        "Matched with DepMap MSIScore (matched): %d" % len(matched),
        "MSI (MSIScore > %.1f): %d" % (msi_threshold, ctr.get("MSI", 0)),
        "MSS: %d" % ctr.get("MSS", 0),
        "",
        "Unmatched / no label:",
    ]
    st = Counter(r["match_status"] for r in mapping_rows)
    for k, v in sorted(st.items()):
        if k != "matched":
            counts_lines.append("  %s: %d" % (k, v))

    counts_path.write_text("\n".join(counts_lines) + "\n")

    # Feasibility narrative
    n_msi, n_mss = ctr.get("MSI", 0), ctr.get("MSS", 0)
    ratio = max(n_msi, n_mss) / max(min(n_msi, n_mss), 1)
    summary_lines = [
        "CCLE MSI vs MSS — mapping summary",
        "=================================",
        "",
        "DepMap files used:",
        "  Model: %s" % model_path.relative_to(root),
        "  Omics signatures: %s" % sig_path.relative_to(root),
        "",
        "Matching logic:",
        "  gene_matrix column name = sample_info 'Cell Line'",
        "  sample_info 'CCLE Code' = DepMap Model 'CCLEName'",
        "  DepMap OmicsGlobalSignatures joined on ModelID (IsDefaultEntryForModel == Yes only)",
        "  MSI if MSIScore > %.1f, else MSS" % msi_threshold,
        "",
        "Counts:",
        "  Proteomics columns (cell lines): %d" % len(matrix_cols),
        "  Successfully labeled (MSI or MSS): %d" % len(matched),
        "  MSI: %d | MSS: %d (ratio max/min ≈ %.2f)" % (n_msi, n_mss, ratio),
        "",
        "Lineage distribution (DepMap OncotreeLineage, matched rows only):",
    ]
    for lin, c in lineage_ctr.most_common():
        summary_lines.append("  %s: %d" % (lin, c))

    summary_lines.extend(
        [
            "",
            "MSI samples by lineage (top):",
        ]
    )
    msi_by_lin = Counter(r["depmap_oncotree_lineage"] for r in matched if r["msi_label"] == "MSI")
    for lin, c in msi_by_lin.most_common(15):
        summary_lines.append("  %s: %d" % (lin, c))

    # Sex distribution (matched rows)
    sex_ct = Counter(((r.get("depmap_sex") or "").strip() or "Unknown") for r in matched)
    summary_lines.extend(
        [
            "",
            "DepMap Sex (matched rows; from Model.csv column Sex):",
        ]
    )
    for s, c in sex_ct.most_common():
        summary_lines.append("  %s: %d" % (s, c))
    summary_lines.append(
        "  Note: Sex is DepMap patient/model annotation, not necessarily assayed in CCLE wet lab."
    )

    # Lineage skew: lineages with only MSI or only MSS
    by_lin = defaultdict(lambda: {"MSI": 0, "MSS": 0})
    for r in matched:
        by_lin[r["depmap_oncotree_lineage"]][r["msi_label"]] += 1
    one_group = [lin for lin, d in by_lin.items() if d["MSI"] == 0 or d["MSS"] == 0]
    summary_lines.extend(
        [
            "",
            "Lineage skew (DepMap OncotreeLineage):",
            "  Lineages with only MSI or only MSS (no within-lineage MSI/MSS mix): %d"
            % len(one_group),
            "  (DA with lineage covariates still identifiable if MSI varies within other lineages;",
            "   see R script for rank checks and rare-lineage collapsing.)",
            "",
            "Feasibility for differential abundance:",
            "  - MSI n=%d vs MSS n=%d — imbalanced; FDR power for MSI effect is limited." % (n_msi, n_mss),
            "  - Pooled across lineages: interpret with **DepMap MSIScore** (not clinical PCR MSI).",
            "  - Prefer **covariate-adjusted** limma (MSI + lineage + sex) per CCLE paper spirit;",
            "    see data/scripts/ccle_DA_msi_vs_mss.R and DA_MSI_vs_MSS/README.txt.",
            "  - Bowel-only subset (~10 MSI vs ~20 MSS here) is a **cleaner** oncologic contrast",
            "    but smaller n — not auto-run unless you add a separate script.",
            "",
            "Unlabeled columns: see match_status in ccle_msi_label_mapping.csv",
        ]
    )

    summary_path = out_dir / "ccle_msi_mapping_summary.txt"
    summary_path.write_text("\n".join(summary_lines) + "\n")

    print("Wrote", map_path.relative_to(root))
    print("Wrote", counts_path.relative_to(root))
    print("Wrote", summary_path.relative_to(root))
    print("Wrote", schema_path.relative_to(root))


if __name__ == "__main__":
    main()
