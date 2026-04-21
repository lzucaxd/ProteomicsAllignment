#!/usr/bin/env bash
# PDC000153 (LUAD / Lung TMT10) — download PSM + run pdc_psm_to_msstatsTMT_protein_matrix.R
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

STUDY="PDC000153"
MANIFEST="${PDC000153_MANIFEST:-$SCRIPT_DIR/manifests/PDC_file_manifest_04112026_120754.csv}"
PSM_ROOT="$SCRIPT_DIR/pdc_psm"
OUT_DIR="$SCRIPT_DIR/results/$STUDY"

# Default .sample.txt: env override > sample_files_msstats_tmt.csv > CPTAC_LOCAL_MIRROR/PDC000153/...
resolve_sample_txt() {
  local raw="${1//\"/}"
  [[ -z "$raw" ]] && return 1
  if [[ -f "$raw" ]]; then echo "$raw"; return 0; fi
  local try="$SCRIPT_DIR/$raw"
  if [[ -f "$try" ]]; then echo "$try"; return 0; fi
  if [[ -n "${CPTAC_LOCAL_MIRROR:-}" ]]; then
    try="${CPTAC_LOCAL_MIRROR%/}/$raw"
    if [[ -f "$try" ]]; then echo "$try"; return 0; fi
  fi
  return 1
}

SAMPLE_TXT="${PDC000153_SAMPLE_TXT:-}"
if [[ -z "$SAMPLE_TXT" && -f "$SCRIPT_DIR/sample_files_msstats_tmt.csv" ]]; then
  RAW=$(awk -F',' -v study="$STUDY" '$1==study {gsub(/^"|"$/,"",$2); print $2; exit}' "$SCRIPT_DIR/sample_files_msstats_tmt.csv")
  if [[ -n "$RAW" ]] && resolved=$(resolve_sample_txt "$RAW"); then
    SAMPLE_TXT="$resolved"
  fi
fi
if [[ -z "$SAMPLE_TXT" && -n "${CPTAC_LOCAL_MIRROR:-}" ]]; then
  try="${CPTAC_LOCAL_MIRROR%/}/PDC000153/CPTAC3_Lung_Adeno_Carcinoma_Proteome.sample.txt"
  [[ -f "$try" ]] && SAMPLE_TXT="$try"
fi

if [[ ! -f "$MANIFEST" ]]; then
  echo "Manifest not found: $MANIFEST"
  exit 1
fi

if [[ ! -f "$SAMPLE_TXT" ]]; then
  echo "Missing CPTAC sample sheet for $STUDY."
  echo "  Set PDC000153_SAMPLE_TXT=/path/to/*.sample.txt, or"
  echo "  export CPTAC_LOCAL_MIRROR=/path/to/parent/of/PDC000153 (see docs/LAB_ONBOARDING.md), or"
  echo "  add a row to data/sample_files_msstats_tmt.csv (path column: PDC000153/your.sample.txt)."
  exit 1
fi

if [[ -x "$SCRIPT_DIR/.venv/bin/python" ]] && "$SCRIPT_DIR/.venv/bin/python" -c "import requests" 2>/dev/null; then
  PYTHON="$SCRIPT_DIR/.venv/bin/python"
else
  PYTHON=python3
fi

N_PSM_IN_MANIFEST="$("$PYTHON" -c "
import csv, sys
n = 0
with open(sys.argv[1], newline='', encoding='utf-8', errors='replace') as f:
    for row in csv.DictReader(f):
        if (row.get('File Name') or '').strip().endswith('.psm'):
            n += 1
print(n)
" "$MANIFEST")"
if [[ "${N_PSM_IN_MANIFEST:-0}" -eq 0 ]]; then
  echo "ERROR: No .psm files in manifest."
  exit 1
fi
echo "[preflight] Manifest lists $N_PSM_IN_MANIFEST .psm file(s)."

echo "[1/2] Download Peptide Spectral Matches (.psm) for $STUDY ..."
if [[ -n "${MAX_PSM_FILES:-}" && "${MAX_PSM_FILES:-0}" != "0" ]]; then
  echo "[info] Limiting download to first $MAX_PSM_FILES PSM files."
  "$PYTHON" "$SCRIPT_DIR/pdc_manifest_downloader.py" \
    --manifest "$MANIFEST" \
    --outdir "$PSM_ROOT" \
    --study-id "$STUDY" \
    --include-category "Peptide Spectral Matches" \
    --ext .psm \
    --max-files "$MAX_PSM_FILES"
else
  "$PYTHON" "$SCRIPT_DIR/pdc_manifest_downloader.py" \
    --manifest "$MANIFEST" \
    --outdir "$PSM_ROOT" \
    --study-id "$STUDY" \
    --include-category "Peptide Spectral Matches" \
    --ext .psm
fi

PSM_DIR="$PSM_ROOT/$STUDY"
if [[ ! -d "$PSM_DIR" ]]; then
  echo "Expected PSM directory not found: $PSM_DIR"
  exit 1
fi

mkdir -p "$OUT_DIR"
echo "[2/2] MSstatsTMT: PSM → protein_summary / gene_matrix ..."
Rscript --no-init-file "$SCRIPT_DIR/pdc_psm_to_msstatsTMT_protein_matrix.R" \
  --psm_dir "$PSM_DIR" \
  --outdir "$OUT_DIR" \
  --sample_txt "$SAMPLE_TXT" \
  --replace_annotation

echo ""
echo "Done. Main outputs:"
echo "  $OUT_DIR/gene_matrix.csv"
echo "  $OUT_DIR/msstats_input.tsv"
echo "  $OUT_DIR/protein_summary.tsv"
