#!/usr/bin/env python
"""
Build subtype-aware mapping for CPTAC breast proteomics study PDC000120.

Inputs (under data/ and results/PDC000120/):
  1) biospecimen/brca_cptac_2020_clinical_data.tsv
  2) biospecimen/S039_BCprospective_01-17_TMT10_Label_to_Sample_Mapping_File_BI_r2.xlsx
  3) results/PDC000120/annotation_filled_corrected.csv
  4) results/PDC000120/gene_matrix.csv

Outputs (all under results/PDC000120/):
  - subtype_mapping_bridge_long.csv
  - subtype_mapping_annotation_to_bridge.csv
  - subtype_mapping_final.csv
  - subtype_mapping_diagnostics.csv
  - subtype_unmatched_samples.csv
  - subtype_ambiguous_matches.csv
  - gene_matrix_subtype_mapping.csv
  - PAM50_tumor_only_samples.csv
  - PAM50_subtype_counts.csv
  - subtype_DA_recommendations.txt

The goal is to produce an analysis-ready mapping from matrix columns to PAM50
and receptor-status labels, with explicit diagnostics for a methods-heavy audience.
"""

import os
import sys
import textwrap
from typing import Tuple, List

import pandas as pd


DATA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES_DIR = os.path.join(DATA_DIR, "results", "PDC000120")
os.makedirs(RES_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# Helper functions: normalization / cleaning
# ---------------------------------------------------------------------------

def norm_colnames(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names: strip, lowercase, replace spaces with underscores."""
    df = df.copy()
    df.columns = (
        df.columns.astype(str)
        .str.replace("\ufeff", "", regex=False)
        .str.strip()
        .str.replace(r"\s+", "_", regex=True)
        .str.lower()
    )
    return df


def norm_channel(x: str) -> str:
    """Normalize TMT channel labels, e.g. '127n' -> '127N', strip spaces."""
    if pd.isna(x):
        return None
    x = str(x).strip()
    # Allow formats like '128C|126' in clinical table: we will treat those later.
    parts = x.split("|")
    normed = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        # Uppercase letter suffix, e.g. 127c -> 127C
        if len(p) >= 3 and p[:3].isdigit():
            head = p[:3]
            tail = p[3:].upper()
            normed.append(head + tail)
        else:
            normed.append(p.upper())
    return "|".join(normed) if normed else None


def norm_plex(x: str) -> str:
    """Normalize TMT plex / mixture names: strip whitespace, keep case-sensitive IDs."""
    if pd.isna(x):
        return None
    return str(x).strip()


def plex_num_from_folder_name(folder_name: str) -> str:
    """
    Extract numeric TMT plex index from CPTAC folder names like:
      '01CPTAC_BCprospective_Proteome_BC_20160911' -> '1'
    Returns None if not parseable.
    """
    if pd.isna(folder_name):
        return None
    s = str(folder_name).strip()
    if len(s) >= 2 and s[:2].isdigit():
        return str(int(s[:2]))  # '01' -> '1'
    return None


def norm_participant(x: str) -> str:
    """Normalize participant / patient IDs, with and without leading X."""
    if pd.isna(x):
        return None
    x = str(x).strip()
    # Many CPTAC IDs are of the form '01BR040' or 'X01BR040'
    if x.startswith("X") and len(x) > 1:
        core = x[1:]
    else:
        core = x
    return core


def add_leading_x(x: str) -> str:
    """Add leading 'X' if not already present, used for matching to clinical table."""
    if pd.isna(x):
        return None
    x = str(x).strip()
    if x.startswith("X"):
        return x
    return "X" + x


def safe_read_tsv(path: str) -> pd.DataFrame:
    """Read TSV with robust header handling."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    return norm_colnames(df)


def safe_read_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    return norm_colnames(df)


def safe_read_gene_matrix(path: str, nrows: int = None) -> pd.DataFrame:
    """Read gene_matrix.csv; it's large so allow restricting rows if needed."""
    df = pd.read_csv(path, dtype=str, nrows=nrows)
    return norm_colnames(df)


# ---------------------------------------------------------------------------
# Step 1: Read all inputs
# ---------------------------------------------------------------------------

def load_inputs() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clinical_path = os.path.join(DATA_DIR, "biospecimen", "brca_cptac_2020_clinical_data.tsv")
    mapping_xlsx = os.path.join(DATA_DIR, "biospecimen", "S039_BCprospective_01-17_TMT10_Label_to_Sample_Mapping_File_BI_r2.xlsx")
    annot_path = os.path.join(RES_DIR, "annotation_filled_corrected.csv")
    matrix_path = os.path.join(RES_DIR, "gene_matrix.csv")

    if not all(os.path.exists(p) for p in [clinical_path, mapping_xlsx, annot_path, matrix_path]):
        missing = [p for p in [clinical_path, mapping_xlsx, annot_path, matrix_path] if not os.path.exists(p)]
        raise FileNotFoundError(f"Missing input files: {missing}")

    print("Reading clinical TSV...")
    clinical = safe_read_tsv(clinical_path)

    print("Reading annotation CSV...")
    annot = safe_read_csv(annot_path)

    print("Reading gene matrix header...")
    matrix = safe_read_gene_matrix(matrix_path, nrows=5)

    print("Reading TMT mapping workbook (all sheets)...")
    # Use engine that supports xlsx
    xls = pd.ExcelFile(mapping_xlsx)
    mapping_sheets = {sheet: norm_colnames(xls.parse(sheet, dtype=str)) for sheet in xls.sheet_names}

    return clinical, annot, matrix, pd.concat(
        [df.assign(_sheet=sheet) for sheet, df in mapping_sheets.items()],
        ignore_index=True
    )


# ---------------------------------------------------------------------------
# Step 2: Unpivot TMT mapping workbook → long bridge
# ---------------------------------------------------------------------------

def unpivot_tmt_mapping(raw_mapping: pd.DataFrame) -> pd.DataFrame:
    """
    Unpivot mapping workbook:
      columns like 'tmt10-126_participant_id', 'tmt10-126_specimen_label', etc.
    Output columns:
      - tmt_plex
      - tmt_channel
      - participant_id
      - specimen_label
      - any useful row-level metadata (e.g. run/plex description)
    """
    df = raw_mapping.copy()

    # This CPTAC workbook is wide, with columns like:
    #   'tmt10-126_participant_id', 'tmt10-126_specimen_label', ...
    # and a per-plex identifier in 'folder_name' (matches proteomics 'Mixture').
    if "_sheet" in df.columns:
        # For PDC000120, our proteomics pipeline uses the proteome mixtures.
        # Keep all sheets for diagnostics, but prioritize Proteome for mapping.
        pass

    if "folder_name" not in df.columns:
        # Some variants may use 'folder' or 'foldername'
        for cand in ("folder", "foldername"):
            if cand in df.columns:
                df["folder_name"] = df[cand]
                break

    wide_cols = [c for c in df.columns if c.startswith("tmt10-") and (c.endswith("_participant_id") or c.endswith("_specimen_label"))]
    if not wide_cols:
        print("WARNING: Could not find expected wide-format TMT10 columns; bridge will be empty.")
        bridge = pd.DataFrame(columns=["tmt_plex", "tmt_channel", "participant_id", "specimen_label", "source_sheet"])
    else:
        rows: List[dict] = []
        # Extract channels present
        channels = sorted({c.split("_")[0].split("-")[-1] for c in wide_cols})
        for _, r in df.iterrows():
            plex = norm_plex(r.get("folder_name"))
            if not plex:
                continue
            plex_num = plex_num_from_folder_name(plex)
            src_sheet = r.get("_sheet")
            for ch in channels:
                pid_col = f"tmt10-{ch}_participant_id"
                spec_col = f"tmt10-{ch}_specimen_label"
                if pid_col not in df.columns and spec_col not in df.columns:
                    continue
                pid = norm_participant(r.get(pid_col)) if pid_col in df.columns else None
                spec = r.get(spec_col) if spec_col in df.columns else None
                if (pid is None or (isinstance(pid, float) and pd.isna(pid))) and (spec is None or (isinstance(spec, float) and pd.isna(spec))):
                    continue
                rows.append({
                    "tmt_plex": plex,
                    "tmt_plex_num": plex_num,
                    "tmt_channel": norm_channel(ch),
                    "participant_id": pid,
                    "specimen_label": str(spec).strip() if pd.notna(spec) else None,
                    "source_sheet": src_sheet,
                    "pcc": r.get("pcc"),
                })
        bridge = pd.DataFrame(rows)
        if not bridge.empty:
            # Prefer Proteome sheet when duplicates exist
            if "source_sheet" in bridge.columns:
                bridge["_sheet_rank"] = bridge["source_sheet"].apply(lambda s: 0 if str(s).lower() == "proteome" else 1)
                bridge = bridge.sort_values(["tmt_plex", "tmt_channel", "_sheet_rank"]).drop_duplicates(["tmt_plex", "tmt_channel"], keep="first")
                bridge = bridge.drop(columns=["_sheet_rank"])
            else:
                bridge = bridge.drop_duplicates(["tmt_plex", "tmt_channel"], keep="first")

    out_path = os.path.join(RES_DIR, "subtype_mapping_bridge_long.csv")
    bridge.to_csv(out_path, index=False)
    print(f"Saved bridge table with {len(bridge)} rows to {out_path}")
    return bridge


# ---------------------------------------------------------------------------
# Step 3: Map annotation file to bridge
# ---------------------------------------------------------------------------

def map_annotation_to_bridge(annot: pd.DataFrame, bridge: pd.DataFrame) -> pd.DataFrame:
    df = annot.copy()
    # Normalize key columns
    for col in ["mixture", "channel", "bioreplicate", "condition"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    df["mixture_norm"] = df.get("mixture", "").apply(norm_plex)
    df["channel_norm"] = df.get("channel", "").apply(norm_channel)

    b = bridge.copy()
    if b.empty:
        print("WARNING: Bridge table is empty; annotation rows will remain unmatched to TMT mapping.")
        # create empty columns so downstream code works
        b["tmt_plex"] = pd.NA
        b["tmt_channel"] = pd.NA
    b["tmt_plex_norm"] = b["tmt_plex"].apply(norm_plex)
    b["tmt_channel_norm"] = b["tmt_channel"].apply(norm_channel)

    # First join on Mixture+Channel
    merged = pd.merge(
        df,
        b,
        left_on=["mixture_norm", "channel_norm"],
        right_on=["tmt_plex_norm", "tmt_channel_norm"],
        how="left",
        suffixes=("_annot", "_bridge"),
    )

    # Diagnostics: matching status
    merged["mapping_status"] = "unmatched"
    # Some mapping workbooks may not provide participant/specimen fields for all entries.
    # Mark as matched if plex+channel matched at all.
    has_plex = merged["tmt_plex"].notna() if "tmt_plex" in merged.columns else False
    has_chan = merged["tmt_channel"].notna() if "tmt_channel" in merged.columns else False
    merged.loc[has_plex & has_chan, "mapping_status"] = "matched_mixture_channel"
    merged["mapping_notes"] = ""

    # Count duplicates or ambiguous mappings (same annot row mapping to >1 bridge row)
    # After our dedupe this shouldn't be frequent, but be explicit:
    # group by Mixture+Channel and see how many bridge entries exist
    bridge_counts = b.groupby(["tmt_plex_norm", "tmt_channel_norm"]).size().reset_index(name="n_bridge_rows")
    merged = pd.merge(
        merged,
        bridge_counts,
        how="left",
        left_on=["mixture_norm", "channel_norm"],
        right_on=["tmt_plex_norm", "tmt_channel_norm"],
        suffixes=("", "_dup"),
    )
    merged["n_bridge_rows"] = merged["n_bridge_rows"].fillna(0).astype(int)
    merged.loc[merged["n_bridge_rows"] > 1, "mapping_status"] = "ambiguous_bridge"
    merged.loc[merged["n_bridge_rows"] > 1, "mapping_notes"] = "Multiple bridge entries for Mixture+Channel"

    out_path = os.path.join(RES_DIR, "subtype_mapping_annotation_to_bridge.csv")
    merged.to_csv(out_path, index=False)
    print(f"Saved annotation→bridge mapping ({len(merged)} rows) to {out_path}")
    return merged


# ---------------------------------------------------------------------------
# Step 4: Map to clinical PAM50 table
# ---------------------------------------------------------------------------

def expand_clinical_multi_channels(clin: pd.DataFrame) -> pd.DataFrame:
    """
    Clinical table sometimes encodes multiple channels / plexes per row, e.g.:
      tmt_channel = '129C|126', tmt_plex = '11|16'
    Expand into one row per (plex, channel).
    """
    df = clin.copy()
    if "tmt_channel" not in df.columns or "tmt_plex" not in df.columns:
        return df

    records = []
    for _, row in df.iterrows():
        chans = str(row["tmt_channel"]).split("|") if pd.notna(row["tmt_channel"]) else [None]
        plexes = str(row["tmt_plex"]).split("|") if pd.notna(row["tmt_plex"]) else [None]
        # align lengths
        if len(plexes) == 1 and len(chans) > 1:
            plexes = plexes * len(chans)
        if len(chans) == 1 and len(plexes) > 1:
            chans = chans * len(plexes)
        for ch, pl in zip(chans, plexes):
            r = row.copy()
            r["tmt_channel"] = norm_channel(ch) if ch is not None else None
            r["tmt_plex"] = norm_plex(pl) if pl is not None else None
            records.append(r)
    return pd.DataFrame(records)


def map_to_clinical(annot_bridge: pd.DataFrame, clinical: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = annot_bridge.copy()

    # Avoid merge suffix collisions: annot_bridge already has bridge columns named tmt_plex/tmt_channel,
    # and we will create key columns tmt_plex_ab/tmt_channel_ab for joining. Rename the bridge columns.
    if "tmt_plex" in df.columns:
        df = df.rename(columns={"tmt_plex": "tmt_plex_bridge"})
    if "tmt_channel" in df.columns:
        df = df.rename(columns={"tmt_channel": "tmt_channel_bridge"})
    if "tmt_plex_num" in df.columns:
        df = df.rename(columns={"tmt_plex_num": "tmt_plex_num_bridge"})

    # Prepare keys in annot+bridge:
    df["tmt_plex_ab"] = df.get("tmt_plex_bridge", "").apply(norm_plex)
    df["tmt_plex_num_ab"] = df.get("tmt_plex_num_bridge", "").apply(lambda x: str(x).strip() if pd.notna(x) and str(x).strip() != "" else None)
    df["tmt_channel_ab"] = df.get("tmt_channel_bridge", "").apply(norm_channel)
    df["participant_core"] = df.get("participant_id", "").apply(norm_participant)

    clin = norm_colnames(clinical)
    # Ensure key columns exist
    for k in ["patientid", "sampleid", "pam50", "tmt_plex", "tmt_channel",
              "er_updated_clinical_status", "pr_clinical_status",
              "erbb2_updated_clinical_status", "tnbc_updated_clinical_status"]:
        if k not in clin.columns:
            # Missing columns will just appear as NaN
            clin[k] = pd.NA

    clin = expand_clinical_multi_channels(clin)
    clin["patient_core"] = clin["patientid"].astype(str).str.strip().apply(norm_participant)
    clin["patientid_with_x"] = clin["patient_core"].apply(add_leading_x)
    clin["tmt_channel"] = clin["tmt_channel"].apply(norm_channel)
    clin["tmt_plex"] = clin["tmt_plex"].apply(norm_plex)
    clin["tmt_plex_num"] = clin["tmt_plex"].astype(str).str.strip()

    # Primary join: Participant + (plex_num, channel)
    merged = pd.merge(
        df,
        clin,
        how="left",
        left_on=["participant_core", "tmt_plex_num_ab", "tmt_channel_ab"],
        right_on=["patient_core", "tmt_plex_num", "tmt_channel"],
        suffixes=("_ab", "_clin"),
    )

    merged["clinical_mapping_status"] = "unmatched"
    matched_mask = merged["pam50"].notna()
    merged.loc[matched_mask, "clinical_mapping_status"] = "matched_participant_plex_channel"

    # Secondary join (still unmatched): by participant only (allows recovery if plex/channel metadata differ)
    unmatched = merged["pam50"].isna()
    if unmatched.any():
        df_unmatched = merged.loc[unmatched].copy()
        clin_pc = clin[["patient_core", "patientid", "sampleid", "pam50",
                        "er_updated_clinical_status", "pr_clinical_status",
                        "erbb2_updated_clinical_status", "tnbc_updated_clinical_status"]].drop_duplicates()
        sec = pd.merge(
            df_unmatched,
            clin_pc,
            how="left",
            left_on=["participant_core"],
            right_on=["patient_core"],
            suffixes=("", "_sec"),
        )
        # Fill PAM50 for cases uniquely identified by plex+channel
        sec_unique = sec.copy()
        # For simplicity, we will fill where exactly one clinical row matches; here we assume dedup above already
        fill_mask = sec_unique["pam50_sec"].notna()
        for col in ["pam50", "er_updated_clinical_status", "pr_clinical_status",
                    "erbb2_updated_clinical_status", "tnbc_updated_clinical_status",
                    "patientid", "sampleid"]:
            src = col + "_sec" if col + "_sec" in sec_unique.columns else col
            if src in sec_unique.columns and col in merged.columns:
                merged.loc[sec_unique.index[fill_mask], col] = sec_unique.loc[fill_mask, src]
        merged.loc[sec_unique.index[fill_mask], "clinical_mapping_status"] = "matched_participant_only"

    # Diagnostics for ambiguity: multiple clinical rows per (plex, channel)
    clin_counts = clin.groupby(["tmt_plex_num", "tmt_channel"]).size().reset_index(name="n_clin_rows")
    merged = pd.merge(
        merged,
        clin_counts,
        how="left",
        left_on=["tmt_plex_num_ab", "tmt_channel_ab"],
        right_on=["tmt_plex_num", "tmt_channel"],
        suffixes=("", "_clin_count"),
    )
    merged["n_clin_rows"] = merged["n_clin_rows"].fillna(0).astype(int)
    merged.loc[(merged["n_clin_rows"] > 1) & merged["pam50"].notna(), "clinical_mapping_status"] = "ambiguous_clinical"

    # Diagnostics tables
    unmatched_rows = merged[merged["pam50"].isna()].copy()
    ambiguous_rows = merged[merged["clinical_mapping_status"].str.contains("ambiguous", na=False)].copy()

    unmatched_out = os.path.join(RES_DIR, "subtype_unmatched_samples.csv")
    ambiguous_out = os.path.join(RES_DIR, "subtype_ambiguous_matches.csv")
    unmatched_rows.to_csv(unmatched_out, index=False)
    ambiguous_rows.to_csv(ambiguous_out, index=False)
    print(f"Saved unmatched samples to {unmatched_out} ({len(unmatched_rows)} rows)")
    print(f"Saved ambiguous matches to {ambiguous_out} ({len(ambiguous_rows)} rows)")

    final_out = os.path.join(RES_DIR, "subtype_mapping_final.csv")
    merged.to_csv(final_out, index=False)
    print(f"Saved final annotation+bridge+clinical mapping to {final_out}")

    diag = merged["clinical_mapping_status"].value_counts(dropna=False).reset_index()
    diag.columns = ["clinical_mapping_status", "n_rows"]
    diag_out = os.path.join(RES_DIR, "subtype_mapping_diagnostics.csv")
    diag.to_csv(diag_out, index=False)
    print(f"Saved clinical mapping diagnostics to {diag_out}")

    return merged, unmatched_rows, ambiguous_rows


# ---------------------------------------------------------------------------
# Step 5: Map to gene matrix columns
# ---------------------------------------------------------------------------

def map_to_gene_matrix(matrix: pd.DataFrame, mapping_final: pd.DataFrame) -> pd.DataFrame:
    mat = matrix.copy()
    # matrix_sample_id columns are all except 'genesymbol' and 'uniprotid'
    non_sample_cols = [c for c in ["genesymbol", "uniprotid"] if c in mat.columns]
    sample_ids = [c for c in mat.columns if c not in non_sample_cols]

    # Annotation's BioReplicate corresponds to matrix column names (matrix is lowercased by norm_colnames)
    mf = mapping_final.copy()
    if "bioreplicate" not in mf.columns:
        print("WARNING: 'bioreplicate' not found in mapping_final; gene-matrix linkage may be incomplete.")
        mf["bioreplicate"] = pd.NA

    br_lower = mf["bioreplicate"].astype(str).str.strip().str.lower()
    sample_ids_set = {str(s).lower() for s in sample_ids}
    lower_to_column = {str(s).lower(): s for s in sample_ids}
    mf["exists_in_gene_matrix"] = br_lower.isin(sample_ids_set)
    # Use actual matrix column name so subsetting works (matrix columns are lowercased)
    mf["matrix_sample_id"] = br_lower.map(lower_to_column).fillna(mf["bioreplicate"])

    # Sample type: from annotation Condition or from clinical if available
    sample_type_col = "condition" if "condition" in mf.columns else None
    sample_types = []
    for _, row in mf.iterrows():
        st = None
        if sample_type_col and pd.notna(row[sample_type_col]):
            st = str(row[sample_type_col]).strip().lower()
        sample_types.append(st)
    mf["sample_type_if_available"] = sample_types

    out = mf.copy()
    out_path = os.path.join(RES_DIR, "gene_matrix_subtype_mapping.csv")
    out.to_csv(out_path, index=False)
    print(f"Saved gene-matrix subtype mapping to {out_path}")
    return out


# ---------------------------------------------------------------------------
# Step 6: Derive tumor-only PAM50 sets
# ---------------------------------------------------------------------------

def derive_tumor_only_sets(gene_matrix_mapping: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = gene_matrix_mapping.copy()

    # Define tumor vs NAT / Norm using sample_type_if_available and clinical ER/PR/HER2 context if needed
    # Here we primarily trust 'condition' (Sample / Norm) at the annotation level and PAM50 presence.
    is_in_matrix = df["exists_in_gene_matrix"].fillna(False)

    # Consider tumor as those with Condition==Sample (case-insensitive) or missing but with PAM50 defined
    cond = df.get("condition", pd.Series([None] * len(df)))
    cond_norm = cond.astype(str).str.lower()
    is_norm = cond_norm.eq("norm")
    is_sample = cond_norm.eq("sample")

    is_tumor = is_sample | ((~is_norm) & df["pam50"].notna())

    tumor_only = df[is_in_matrix & is_tumor & df["pam50"].notna()].copy()

    # Basic subtype counts
    pam50_counts = tumor_only["pam50"].value_counts(dropna=False).reset_index()
    pam50_counts.columns = ["pam50", "n_tumor_samples"]

    tumor_out = os.path.join(RES_DIR, "PAM50_tumor_only_samples.csv")
    counts_out = os.path.join(RES_DIR, "PAM50_subtype_counts.csv")
    tumor_only.to_csv(tumor_out, index=False)
    pam50_counts.to_csv(counts_out, index=False)
    print(f"Saved tumor-only PAM50 sample table to {tumor_out} ({len(tumor_only)} rows)")
    print(f"Saved PAM50 subtype counts to {counts_out}")

    return tumor_only, pam50_counts


# ---------------------------------------------------------------------------
# Step 6b: DA-ready sample table (one row per gene_matrix column)
# ---------------------------------------------------------------------------

def build_da_sample_annotation(
    gene_matrix_mapping: pd.DataFrame,
    matrix: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build tables for DA: one row per gene_matrix column.
    - DA_sample_annotation: all matrix columns with sample_type and PAM50 when available.
    - DA_subtype_tumor_only: tumor columns with non-missing PAM50 (for subtype DA).
    Uses biospecimen for sample_type (Primary Tumor vs Solid Tissue Normal).
    """
    non_sample_cols = [c for c in ["genesymbol", "uniprotid"] if c in matrix.columns]
    sample_ids = [c for c in matrix.columns if c not in non_sample_cols]

    # Biospecimen: aliquot -> Sample Type
    bio_path = os.path.join(DATA_DIR, "biospecimen", "PDC_study_biospecimen_03162026_190026.csv")
    aliquot_to_type = {}
    if os.path.exists(bio_path):
        bio = pd.read_csv(bio_path, dtype=str)
        bio.columns = [c.replace("\ufeff", "").strip() for c in bio.columns]
        aid_col = next((c for c in bio.columns if "aliquot" in c.lower() and "submitter" in c.lower()), None)
        st_col = next((c for c in bio.columns if "sample" in c.lower() and "type" in c.lower()), None)
        if aid_col and st_col:
            for _, r in bio.iterrows():
                aid = str(r[aid_col]).strip()
                st = str(r[st_col]).strip()
                if aid:
                    # store lower-cased key so lookups can be case-insensitive
                    aliquot_to_type[aid.lower()] = st

    # One row per unique matrix column; prefer rows that have PAM50 when deduplicating
    key_cols = ["bioreplicate", "matrix_sample_id", "pam50", "patientid", "sampleid",
                "er_updated_clinical_status", "pr_clinical_status",
                "erbb2_updated_clinical_status", "tnbc_updated_clinical_status",
                "condition", "mixture"]
    have = [c for c in key_cols if c in gene_matrix_mapping.columns]
    sub = gene_matrix_mapping[have + ["exists_in_gene_matrix"]].copy()
    sub = sub[sub["exists_in_gene_matrix"].fillna(False)]
    sub = sub.sort_values(by="pam50", na_position="last")  # rows with PAM50 first
    sub = sub.drop_duplicates(subset=["matrix_sample_id"], keep="first")

    # Map sample_type from biospecimen
    def sample_type_for_da(row):
        aid = row.get("bioreplicate") or row.get("matrix_sample_id")
        if pd.isna(aid):
            return "Unknown"
        # case-insensitive lookup against biospecimen Aliquot Submitter ID
        st = aliquot_to_type.get(str(aid).strip().lower())
        if pd.isna(st) or st == "":
            cond = str(row.get("condition", "")).strip().lower()
            if cond == "norm":
                return "Norm"
            return "Unknown"
        if "primary tumor" in str(st).lower():
            return "Tumor"
        if "solid tissue normal" in str(st).lower():
            return "NAT"
        return st

    sub["sample_type"] = sub.apply(sample_type_for_da, axis=1)
    # Keep matrix_sample_id as actual matrix column name (for subsetting and missing check)

    # Full annotation: every matrix column that we could map
    da_annot = sub.copy()
    # Ensure all matrix columns appear (even if not in mapping)
    missing = set(sample_ids) - set(da_annot["matrix_sample_id"])
    if missing:
        extra = pd.DataFrame({
            "matrix_sample_id": list(missing),
            "bioreplicate": list(missing),
            "sample_type": [aliquot_to_type.get(str(a).strip().lower(), "Unknown") for a in missing],
            "pam50": pd.NA,
            "patientid": pd.NA,
            "exists_in_gene_matrix": True,
        })
        for c in have:
            if c not in extra.columns:
                extra[c] = pd.NA
        da_annot = pd.concat([da_annot, extra], ignore_index=True)

    # Subtype tumor-only: Tumor + non-missing PAM50 (for Basal vs LumA etc.)
    da_tumor = da_annot[
        (da_annot["sample_type"] == "Tumor") &
        da_annot["pam50"].notna() &
        (da_annot["pam50"].astype(str).str.strip() != "")
    ].copy()

    # Subtype counts (unique matrix columns)
    subtype_counts = da_tumor["pam50"].value_counts(dropna=False).reset_index()
    subtype_counts.columns = ["PAM50", "n"]

    out_annot = os.path.join(RES_DIR, "DA_sample_annotation.csv")
    out_tumor = os.path.join(RES_DIR, "DA_subtype_tumor_only.csv")
    out_counts = os.path.join(RES_DIR, "DA_subtype_counts.csv")

    cols_annot = ["matrix_sample_id", "sample_type", "pam50", "patientid", "sampleid",
                  "er_updated_clinical_status", "pr_clinical_status",
                  "erbb2_updated_clinical_status", "tnbc_updated_clinical_status", "mixture"]
    cols_annot = [c for c in cols_annot if c in da_annot.columns]
    da_annot[cols_annot].to_csv(out_annot, index=False)
    da_tumor[cols_annot].to_csv(out_tumor, index=False)
    subtype_counts.to_csv(out_counts, index=False)

    print(f"  DA_sample_annotation.csv: {len(da_annot)} rows (one per matrix column with mapping)")
    print(f"  DA_subtype_tumor_only.csv: {len(da_tumor)} rows (Tumor + PAM50, for subtype DA)")
    print(f"  DA_subtype_counts.csv: counts by PAM50 for subtype DA")

    return da_annot, da_tumor, subtype_counts


# ---------------------------------------------------------------------------
# Step 7: Suggest subtype contrasts
# ---------------------------------------------------------------------------

def suggest_contrasts(pam50_counts: pd.DataFrame, min_per_group: int = 5) -> List[str]:
    counts = pam50_counts.set_index("pam50")["n_tumor_samples"].to_dict()
    # Canonical subtypes
    n_basal = counts.get("Basal", 0) + counts.get("Basal-like", 0)
    n_luma = counts.get("LumA", 0)
    n_lumb = counts.get("LumB", 0)
    n_her2 = counts.get("Her2", 0) + counts.get("HER2", 0)
    n_norm_like = counts.get("Normal-like", 0)

    recs = []

    def add_if_ok(name: str, n1: int, n2: int):
        if n1 >= min_per_group and n2 >= min_per_group:
            recs.append(f"{name}: {n1} vs {n2} samples (both ≥ {min_per_group})")

    add_if_ok("Basal vs LumA", n_basal, n_luma)
    add_if_ok("Basal vs Luminal (LumA+LumB)", n_basal, n_luma + n_lumb)
    add_if_ok("Her2 vs LumA", n_her2, n_luma)
    add_if_ok("LumA vs LumB", n_luma, n_lumb)
    if n_norm_like >= min_per_group and n_luma >= min_per_group:
        add_if_ok("Normal-like vs LumA", n_norm_like, n_luma)

    if not recs:
        recs.append("No subtype contrasts meet the minimum per-group size threshold.")
    return recs


def write_recommendations(pam50_counts: pd.DataFrame, recommendations: List[str]) -> None:
    out_path = os.path.join(RES_DIR, "subtype_DA_recommendations.txt")
    lines = []
    lines.append("Subtype DA recommendations for CPTAC BRCA (PDC000120)")
    lines.append("")
    lines.append("PAM50 subtype counts (tumor-only, matrix-present):")
    for _, row in pam50_counts.iterrows():
        lines.append(f"  - {row['pam50']}: {row['n_tumor_samples']}")
    lines.append("")
    lines.append("Recommended contrasts (based on counts):")
    for r in recommendations:
        lines.append(f"  - {r}")
    lines.append("")
    lines.append("DA-ready tables (one row per gene_matrix column):")
    lines.append("  - DA_sample_annotation.csv: all matrix columns with sample_type (Tumor/NAT) and PAM50 when available.")
    lines.append("  - DA_subtype_tumor_only.csv: tumor columns with PAM50 only; use this for subtype DA (e.g. Basal vs LumA).")
    lines.append("  - DA_subtype_counts.csv: counts by PAM50 for tumor-only samples.")
    lines.append("  Run: Rscript scripts/DA_subtype_CPTAC_breast.R (from data/) for subtype limma and volcano plots.")
    lines.append("")
    lines.append("Caveats:")
    lines.append("  - These recommendations assume PAM50 calls in the clinical file are trustworthy.")
    lines.append("  - Some subtypes (e.g., Her2, Normal-like) may be underpowered depending on thresholds.")
    lines.append("  - For methods-heavy use, consider modeling subtype as a multi-level factor with limma or MSstatsTMT and extracting specific contrasts.")
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"Saved subtype DA recommendations to {out_path}")


# ---------------------------------------------------------------------------
# Step 8: Final console summary
# ---------------------------------------------------------------------------

def print_final_summary(annot: pd.DataFrame,
                        bridge: pd.DataFrame,
                        mapping_final: pd.DataFrame,
                        gene_matrix_mapping: pd.DataFrame,
                        tumor_only: pd.DataFrame) -> None:
    n_annot = len(annot)
    n_bridge = len(bridge)
    n_mapped_bridge = (mapping_final["mapping_status"] != "unmatched").sum() if "mapping_status" in mapping_final.columns else pd.NA
    n_with_pam50 = mapping_final["pam50"].notna().sum() if "pam50" in mapping_final.columns else pd.NA
    n_in_matrix = gene_matrix_mapping["exists_in_gene_matrix"].sum()
    n_tumor_only = len(tumor_only)

    summary = textwrap.dedent(f"""
    === Subtype mapping summary (PDC000120) ===
      Annotation rows:                     {n_annot}
      Unique TMT plex/channel in bridge:   {n_bridge}
      Annotation rows matched to bridge:   {n_mapped_bridge}
      Annotation rows with PAM50 label:    {n_with_pam50}
      Samples present in gene_matrix:      {n_in_matrix}
      Tumor-only PAM50 samples (matrix):   {n_tumor_only}
    """).strip()
    print(summary)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    clinical, annot, matrix, mapping_raw = load_inputs()

    bridge = unpivot_tmt_mapping(mapping_raw)
    annot_bridge = map_annotation_to_bridge(annot, bridge)
    mapping_final, unmatched_rows, ambiguous_rows = map_to_clinical(annot_bridge, clinical)
    gene_matrix_mapping = map_to_gene_matrix(matrix, mapping_final)
    tumor_only, pam50_counts = derive_tumor_only_sets(gene_matrix_mapping)
    da_annot, da_tumor, da_subtype_counts = build_da_sample_annotation(gene_matrix_mapping, matrix)
    recs = suggest_contrasts(da_subtype_counts.rename(columns={"PAM50": "pam50", "n": "n_tumor_samples"}))
    write_recommendations(da_subtype_counts.rename(columns={"PAM50": "pam50", "n": "n_tumor_samples"}), recs)
    print_final_summary(annot, bridge, mapping_final, gene_matrix_mapping, tumor_only)


if __name__ == "__main__":
    sys.exit(main())

