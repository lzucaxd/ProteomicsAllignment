#!/usr/bin/env bash
# =============================================================================
# Batch runner: download PSMs, run MSstatsTMT, keep only final outputs.
#
# Space-efficient: processes one study at a time, deletes intermediates
# (parsed_psm_long.tsv, msstats_input.tsv, protein_summary.tsv) and raw
# PSM downloads after each study finishes successfully.
#
# Usage:
#   cd data/ && ./run_batch_studies.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Directory that contains PDC000153/, PDC000204/, ... (each with *.sample.txt).
# Same layout as a local CPTAC "data" tree. Required for the paths below.
if [[ -z "${CPTAC_LOCAL_MIRROR:-}" ]]; then
  echo "ERROR: export CPTAC_LOCAL_MIRROR=/path/to/dir_containing_PDC_study_folders"
  echo "See docs/LAB_ONBOARDING.md"
  exit 1
fi
MIRROR="${CPTAC_LOCAL_MIRROR%/}"

PSM_ROOT="$SCRIPT_DIR/pdc_psm"
INTERMEDIATES="parsed_psm_long.tsv msstats_input.tsv protein_summary.tsv"

if [[ -x "$SCRIPT_DIR/.venv/bin/python" ]] && "$SCRIPT_DIR/.venv/bin/python" -c "import requests" 2>/dev/null; then
  PYTHON="$SCRIPT_DIR/.venv/bin/python"
else
  PYTHON=python3
fi

run_one_study() {
  local STUDY="$1"
  local MANIFEST_REL="$2"
  local SAMPLE_TXT="$3"

  local MANIFEST="$SCRIPT_DIR/$MANIFEST_REL"
  local OUT_DIR="$SCRIPT_DIR/results/$STUDY"

  echo ""
  echo "================================================================"
  echo "  $STUDY"
  echo "================================================================"

  if [[ -f "$OUT_DIR/gene_matrix.csv" ]]; then
    echo "[skip] $STUDY — gene_matrix.csv already exists in $OUT_DIR"
    return 0
  fi

  if [[ ! -f "$MANIFEST" ]]; then
    echo "[ERROR] Manifest not found: $MANIFEST — skipping $STUDY"
    return 0
  fi
  if [[ ! -f "$SAMPLE_TXT" ]]; then
    echo "[ERROR] sample.txt not found: $SAMPLE_TXT — skipping $STUDY"
    return 0
  fi

  local N_PSM
  N_PSM=$("$PYTHON" -c "
import csv, sys
n = 0
with open(sys.argv[1], newline='', encoding='utf-8', errors='replace') as f:
    r = csv.DictReader(f)
    nc = 'File Name' if 'File Name' in (r.fieldnames or []) else None
    if not nc: print(0); sys.exit(0)
    for row in r:
        if (row.get(nc) or '').strip().endswith('.psm'): n += 1
print(n)
" "$MANIFEST")
  if [[ "${N_PSM:-0}" -eq 0 ]]; then
    echo "[ERROR] Manifest has 0 .psm files — skipping $STUDY"
    return 0
  fi
  echo "[preflight] $N_PSM .psm files in manifest."

  echo "[1/3] Downloading PSMs for $STUDY ..."
  "$PYTHON" "$SCRIPT_DIR/pdc_manifest_downloader.py" \
    --manifest "$MANIFEST" \
    --outdir "$PSM_ROOT" \
    --study-id "$STUDY" \
    --include-category "Peptide Spectral Matches" \
    --ext .psm

  local PSM_DIR="$PSM_ROOT/$STUDY"
  if [[ ! -d "$PSM_DIR" ]]; then
    echo "[ERROR] PSM directory not found after download: $PSM_DIR — skipping $STUDY"
    return 0
  fi

  mkdir -p "$OUT_DIR"
  echo "[2/3] Running MSstatsTMT pipeline ..."
  if ! Rscript --no-init-file "$SCRIPT_DIR/pdc_psm_to_msstatsTMT_protein_matrix.R" \
    --psm_dir "$PSM_DIR" \
    --outdir "$OUT_DIR" \
    --sample_txt "$SAMPLE_TXT" \
    --replace_annotation; then
    echo "[ERROR] MSstatsTMT pipeline failed for $STUDY — see $OUT_DIR for partial outputs"
    return 0
  fi

  if [[ ! -f "$OUT_DIR/gene_matrix.csv" ]]; then
    echo "[ERROR] gene_matrix.csv not produced for $STUDY — keeping intermediates for debugging"
    return 0
  fi

  echo "[3/3] Cleaning up intermediates for $STUDY ..."
  for f in $INTERMEDIATES; do
    rm -f "$OUT_DIR/$f"
  done
  rm -rf "$PSM_DIR"
  echo "[done] $STUDY — kept gene_matrix.csv, protein_matrix_wide.csv, annotations, plots, qc_summary"
}

# ---------------------------------------------------------------------------
# Study list: STUDY_ID  MANIFEST  SAMPLE_TXT
# ---------------------------------------------------------------------------
run_one_study PDC000153 \
  "manifests/PDC_file_manifest_04112026_120754.csv" \
  "$MIRROR/PDC000153/CPTAC3_Lung_Adeno_Carcinoma_Proteome.sample.txt"

run_one_study PDC000204 \
  "manifests/PDC_file_manifest_04112026_120906.csv" \
  "$MIRROR/PDC000204/CPTAC3_Glioblastoma_Multiforme_Proteome.sample.txt"

run_one_study PDC000221 \
  "manifests/PDC_file_manifest_04112026_120953.csv" \
  "$MIRROR/PDC000221/CPTAC3_Head_and_Neck_Carcinoma_Proteome.sample.txt"

run_one_study PDC000234 \
  "manifests/PDC_file_manifest_04112026_121016.csv" \
  "$MIRROR/PDC000234/CPTAC3_Lung_Squamous_Cell_Carcinoma_Proteome.sample.txt"

run_one_study PDC000270 \
  "manifests/PDC_file_manifest_04112026_121046.csv" \
  "$MIRROR/PDC000270/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.sample.txt"

run_one_study PDC000464 \
  "manifests/PDC_file_manifest_04112026_121119.csv" \
  "$MIRROR/PDC000464/CPTAC3_non-ccRCC_JHU_Proteome.sample.txt"

echo ""
echo "========================================"
echo "  Batch complete. Final outputs:"
echo "========================================"
for STUDY in PDC000153 PDC000204 PDC000221 PDC000234 PDC000270 PDC000464; do
  OUT="$SCRIPT_DIR/results/$STUDY"
  if [[ -f "$OUT/gene_matrix.csv" ]]; then
    rows=$(wc -l < "$OUT/gene_matrix.csv" | tr -d ' ')
    echo "  $STUDY: gene_matrix.csv ($rows rows)"
  else
    echo "  $STUDY: MISSING"
  fi
done
