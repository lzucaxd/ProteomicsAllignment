"""Build standardized sample and feature metadata for benchmark tasks."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _norm_cell_line_token(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def _ccle_column_by_line_norm(ccle_columns: list[str]) -> dict[str, str]:
    """Map normalized line token -> first matching CCLE matrix column name."""
    out: dict[str, str] = {}
    for col in ccle_columns:
        key = _norm_cell_line_token(col)
        out.setdefault(key, col)
    return out


def average_cal120_ccle_columns(ccle: pd.DataFrame, out_name: str = "CAL120_BREAST") -> pd.DataFrame:
    """
    Collapse every column whose normalized name is CAL120 into one column (row-wise mean).
    If only one such column exists, rename it to out_name. Required so CAL120 is one
    biological unit (two plexes in the annotation file, one column in analysis).
    """
    df = ccle.copy()
    cols = list(df.columns)
    cal_cols = [c for c in cols if _norm_cell_line_token(c) == "CAL120"]
    if len(cal_cols) >= 2:
        df[out_name] = df[cal_cols].mean(axis=1, skipna=True)
        df = df.drop(columns=cal_cols)
        logger.info("CCLE: averaged %d CAL-120 columns into %s", len(cal_cols), out_name)
    elif len(cal_cols) == 1:
        df = df.rename(columns={cal_cols[0]: out_name})
        logger.info("CCLE: renamed single CAL-120 column %s -> %s", cal_cols[0], out_name)
    return df


def _resolve_ccle_column_for_cell_line(
    cell_line: str,
    ccle_columns: list[str],
    col_lookup: dict[str, str],
    cal120_breast_col: str = "CAL120_BREAST",
) -> str | None:
    """Map annotation cell_line to a column in the (possibly CAL120-merged) CCLE matrix."""
    if _norm_cell_line_token(cell_line) == "CAL120":
        if cal120_breast_col in ccle_columns:
            return cal120_breast_col
    key = _norm_cell_line_token(cell_line)
    return col_lookup.get(key)


def _ccle_rows_from_subtype_annotation_csv(
    ccle: pd.DataFrame,
    csv_path: Path,
    group_col: str = "BvL_group",
    processed_path: Path | None = None,
) -> list[dict[str, Any]]:
    """
    Load CCLE breast subtype annotations. Only Basal and Luminal (HER2 excluded).
    Prefers processed one-row-per-cell_line table if present.
    """
    path_use = processed_path if processed_path is not None and processed_path.is_file() else csv_path
    ann = pd.read_csv(path_use)
    if group_col not in ann.columns:
        raise ValueError(f"Annotation CSV missing column {group_col!r}: {path_use}")
    if "cell_line" not in ann.columns:
        raise ValueError(f"Annotation CSV missing 'cell_line': {path_use}")

    ann = ann[ann[group_col].isin(["Basal", "Luminal"])].copy()
    # processed file is already deduped; raw v2 needs dedup
    if "n_plexes" not in ann.columns:
        ann = ann.drop_duplicates(subset=["cell_line"], keep="first")

    ccols = list(ccle.columns)
    col_lookup = _ccle_column_by_line_norm(ccols)
    rows: list[dict[str, Any]] = []
    for _, row in ann.iterrows():
        line = str(row["cell_line"])
        grp = str(row[group_col])
        condition = "Basal" if grp == "Basal" else "Luminal"
        sid = _resolve_ccle_column_for_cell_line(line, ccols, col_lookup)
        if sid is None:
            logger.warning("Subtype annotation: no CCLE column for cell_line=%s", line)
            continue
        rows.append({"sample_id": sid, "domain": "CCLE", "condition": condition, "study_id": "CCLE"})
    return rows


def _ccle_breast_rows_from_full_v2_annotation(
    ccle: pd.DataFrame,
    csv_path: Path,
) -> list[dict[str, Any]]:
    """All breast lines in v2 (Basal + Luminal + HER2) for BvL CCLE Breast arm."""
    ann = pd.read_csv(csv_path)
    if "cell_line" not in ann.columns or "BvL_group" not in ann.columns:
        raise ValueError(f"Expected cell_line and BvL_group in {csv_path}")
    ann = ann.drop_duplicates(subset=["cell_line"], keep="first")
    ccols = list(ccle.columns)
    col_lookup = _ccle_column_by_line_norm(ccols)
    rows: list[dict[str, Any]] = []
    for _, row in ann.iterrows():
        line = str(row["cell_line"])
        sid = _resolve_ccle_column_for_cell_line(line, ccols, col_lookup)
        if sid is None:
            logger.warning("BvL breast annotation: no CCLE column for cell_line=%s", line)
            continue
        rows.append({"sample_id": sid, "domain": "CCLE", "condition": "Breast", "study_id": "CCLE"})
    return rows


def build_sample_meta(
    cptac_matrices: dict[str, pd.DataFrame],
    ccle_matrix: pd.DataFrame,
    task_config: dict[str, Any],
    ccle_sample_info: pd.DataFrame | None = None,
    repo_root: Path | None = None,
) -> pd.DataFrame:
    """
    Build a standardized sample metadata DataFrame for a benchmark task.

    Columns: sample_id, domain, condition, study_id
    """
    task_name = task_config["task_name"]

    if task_name == "breast_subtype":
        return _build_subtype_meta(cptac_matrices, ccle_matrix, task_config, repo_root=repo_root)
    elif task_name == "breast_vs_lung":
        return _build_bvl_meta(
            cptac_matrices, ccle_matrix, task_config, ccle_sample_info, repo_root=repo_root
        )
    else:
        raise ValueError(f"Unknown task: {task_name}")


def _build_subtype_meta(
    cptac: dict[str, pd.DataFrame],
    ccle: pd.DataFrame,
    cfg: dict[str, Any],
    repo_root: Path | None = None,
) -> pd.DataFrame:
    """Breast subtype: Luminal vs Basal from CPTAC and CCLE."""
    rows = []

    # CPTAC samples from subtype mapping
    mapping_path = cfg["cptac"].get("subtype_mapping")
    if mapping_path:
        sm = pd.read_csv(mapping_path)
        st_col = "sample_type" if "sample_type" in sm.columns else "sample_type_if_available"
        study_id = cfg["cptac"]["studies"][0]
        if study_id in cptac:
            mat_cols = list(cptac[study_id].columns)
            mat_cols_lower = {c.lower(): c for c in mat_cols}

            luminal_pam50 = [v.lower() for v in cfg["cptac"]["conditions"].get("Luminal", [])]
            basal_pam50 = [v.lower() for v in cfg["cptac"]["conditions"].get("Basal", [])]

            for _, row in sm.iterrows():
                if str(row.get(st_col, "")).lower() != "sample":
                    continue
                if row.get("exists_in_gene_matrix") is not True and row.get("exists_in_gene_matrix") != "TRUE":
                    continue
                pam50 = str(row.get("pam50", "")).lower()
                if pam50 in luminal_pam50:
                    condition = "Luminal"
                elif pam50 in basal_pam50:
                    condition = "Basal"
                else:
                    continue
                sid = str(row.get("matrix_sample_id", ""))
                matched = mat_cols_lower.get(sid.lower())
                if matched:
                    rows.append({"sample_id": matched, "domain": "CPTAC", "condition": condition, "study_id": study_id})

    # CCLE samples — prefer curated annotation CSV (v2) for power and consistency
    ccle_cfg = cfg.get("ccle", {})
    ann_rel = ccle_cfg.get("subtype_annotation_csv")
    ann_path: Path | None = None
    if ann_rel and repo_root is not None:
        ann_path = Path(ann_rel)
        if not ann_path.is_absolute():
            ann_path = repo_root / ann_path

    if ann_path is not None and ann_path.is_file():
        group_col = ccle_cfg.get("subtype_group_column", "BvL_group")
        proc_rel = ccle_cfg.get("subtype_annotation_processed_csv", "data/processed/ccle_breast_subtype_annotation_processed.csv")
        proc_path = Path(proc_rel)
        if repo_root is not None and not proc_path.is_absolute():
            proc_path = repo_root / proc_path
        processed_arg = proc_path if proc_path.is_file() else None
        n_before = len(rows)
        rows.extend(
            _ccle_rows_from_subtype_annotation_csv(
                ccle, ann_path, group_col=group_col, processed_path=processed_arg
            )
        )
        logger.info(
            "CCLE subtype: added %d samples from annotation %s",
            len(rows) - n_before,
            ann_path,
        )
    else:
        if ann_rel:
            logger.warning("subtype_annotation_csv not found (%s); using basal_lines/luminal_lines", ann_path)
        ccle_cols = list(ccle.columns)
        for line_name in ccle_cfg.get("basal_lines", []):
            pat = re.escape(line_name).replace(r"\-", ".")
            matches = [c for c in ccle_cols if re.search(pat, c, re.IGNORECASE)]
            if matches:
                rows.append({"sample_id": matches[0], "domain": "CCLE", "condition": "Basal", "study_id": "CCLE"})
        for line_name in ccle_cfg.get("luminal_lines", []):
            pat = re.escape(line_name).replace(r"\-", ".")
            matches = [c for c in ccle_cols if re.search(pat, c, re.IGNORECASE)]
            if matches:
                rows.append({"sample_id": matches[0], "domain": "CCLE", "condition": "Luminal", "study_id": "CCLE"})

    return pd.DataFrame(rows)


def _build_bvl_meta(
    cptac: dict[str, pd.DataFrame],
    ccle: pd.DataFrame,
    cfg: dict[str, Any],
    ccle_info: pd.DataFrame | None = None,
    repo_root: Path | None = None,
) -> pd.DataFrame:
    """Breast vs Lung from CPTAC studies + CCLE (breast lines from v2 optional, lung from tissue map)."""
    rows = []

    breast_study = cfg["cptac"].get("breast_study", "PDC000120")
    lung_study = cfg["cptac"].get("lung_study", "PDC000153")

    if breast_study in cptac:
        for sid in cptac[breast_study].columns:
            rows.append({"sample_id": sid, "domain": "CPTAC", "condition": "Breast", "study_id": breast_study})
    if lung_study in cptac:
        for sid in cptac[lung_study].columns:
            rows.append({"sample_id": sid, "domain": "CPTAC", "condition": "Lung", "study_id": lung_study})

    ccle_cfg = cfg.get("ccle", {})
    breast_ann_rel = ccle_cfg.get("breast_lines_annotation_csv")
    breast_ann_path: Path | None = None
    if breast_ann_rel and repo_root is not None:
        breast_ann_path = Path(breast_ann_rel)
        if not breast_ann_path.is_absolute():
            breast_ann_path = repo_root / breast_ann_path

    use_ann_breast = breast_ann_path is not None and breast_ann_path.is_file()

    if use_ann_breast:
        rows.extend(_ccle_breast_rows_from_full_v2_annotation(ccle, breast_ann_path))

    # CCLE lung (and breast if not using annotation) via tissue mapping
    if ccle_info is not None and "Tissue of Origin" in ccle_info.columns:
        tissue_map = dict(zip(ccle_info["Cell Line"], ccle_info["Tissue of Origin"]))
        normalize = lambda x: re.sub(r"[^A-Za-z0-9]", "", x).upper()
        norm_map = {normalize(k): (k, v) for k, v in tissue_map.items()}

        for col in ccle.columns:
            tissue = None
            if col in tissue_map:
                tissue = tissue_map[col]
            else:
                norm_col = normalize(col)
                if norm_col in norm_map:
                    tissue = norm_map[norm_col][1]
            if not tissue:
                continue
            tl = tissue.lower()
            if tl == "lung":
                rows.append({
                    "sample_id": col,
                    "domain": "CCLE",
                    "condition": "Lung",
                    "study_id": "CCLE",
                })
            elif tl == "breast" and not use_ann_breast:
                rows.append({
                    "sample_id": col,
                    "domain": "CCLE",
                    "condition": "Breast",
                    "study_id": "CCLE",
                })

    return pd.DataFrame(rows)


def build_feature_meta(
    shared_genes: list[str],
    all_genes: list[str] | None = None,
) -> pd.DataFrame:
    """Build feature metadata with inclusion/exclusion status."""
    if all_genes is None:
        all_genes = shared_genes

    records = []
    shared_set = set(shared_genes)
    for g in all_genes:
        records.append({
            "gene": g,
            "included": g in shared_set,
            "exclusion_reason": "" if g in shared_set else "not_in_intersection",
        })
    return pd.DataFrame(records)
