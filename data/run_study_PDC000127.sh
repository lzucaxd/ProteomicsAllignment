#!/usr/bin/env bash
# PDC000127 (ccRCC TMT10) — download PSM + run pdc_psm_to_msstatsTMT_protein_matrix.R
# Run from repo root OR from data/ (script cd's to data/).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

STUDY="PDC000127"
MANIFEST="${PDC000127_MANIFEST:-$SCRIPT_DIR/manifests/PDC_file_manifest_04072026_160622.csv}"
PSM_ROOT="$SCRIPT_DIR/pdc_psm"
OUT_DIR="$SCRIPT_DIR/results/$STUDY"
DEFAULT_SAMPLE="$SCRIPT_DIR/cptac_samples/$STUDY/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.sample.txt"
SAMPLE_TXT="${PDC000127_SAMPLE_TXT:-$DEFAULT_SAMPLE}"

if [[ ! -f "$MANIFEST" ]]; then
  echo "Manifest not found: $MANIFEST"
  exit 1
fi

if [[ ! -f "$SAMPLE_TXT" ]]; then
  echo "Missing CPTAC sample sheet for MSstatsTMT annotation:"
  echo "  $SAMPLE_TXT"
  echo ""
  echo "Copy or symlink CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.sample.txt here."
  echo "See: $SCRIPT_DIR/cptac_samples/$STUDY/README.md"
  exit 1
fi

if [[ -x "$SCRIPT_DIR/.venv/bin/python" ]] && "$SCRIPT_DIR/.venv/bin/python" -c "import requests" 2>/dev/null; then
  PYTHON="$SCRIPT_DIR/.venv/bin/python"
else
  PYTHON=python3
fi

# Pipeline parses tab-delimited .psm only (not mzIdentML .mzid.gz).
N_PSM_IN_MANIFEST="$("$PYTHON" -c "
import csv, sys
path = sys.argv[1]
n = 0
with open(path, newline='', encoding='utf-8', errors='replace') as f:
    r = csv.DictReader(f)
    name_col = 'File Name' if 'File Name' in (r.fieldnames or []) else None
    if not name_col:
        print(0)
        sys.exit(0)
    for row in r:
        fn = (row.get(name_col) or '').strip()
        if fn.endswith('.psm'):
            n += 1
print(n)
" "$MANIFEST")"
if [[ "${N_PSM_IN_MANIFEST:-0}" -eq 0 ]]; then
  echo "ERROR: This manifest has no .psm (tabular PSM) files — found 0 filenames ending in .psm."
  echo "Your export is likely mzIdentML only (.mzid.gz). This repo pipeline needs CPTAC Text PSM files."
  echo ""
  echo "In the PDC portal: Files → Peptide Spectral Matches → filter File Type = **Text** (or files ending in .psm),"
  echo "export a new manifest CSV, then:"
  echo "  PDC000127_MANIFEST=$SCRIPT_DIR/manifests/<your_new_manifest>.csv ./run_study_PDC000127.sh"
  echo "See: $SCRIPT_DIR/cptac_samples/$STUDY/README.md"
  exit 1
fi
echo "[preflight] Manifest lists $N_PSM_IN_MANIFEST .psm file(s)."

if [[ -n "${MAX_PSM_FILES:-}" && "$MAX_PSM_FILES" != "0" ]]; then
  echo "[info] Limiting download to first $MAX_PSM_FILES PSM files (MAX_PSM_FILES)."
  echo "[1/2] Download Peptide Spectral Matches (.psm) for $STUDY ..."
  "$PYTHON" "$SCRIPT_DIR/pdc_manifest_downloader.py" \
    --manifest "$MANIFEST" \
    --outdir "$PSM_ROOT" \
    --study-id "$STUDY" \
    --include-category "Peptide Spectral Matches" \
    --ext .psm \
    --max-files "$MAX_PSM_FILES"
else
  echo "[1/2] Download Peptide Spectral Matches (.psm) for $STUDY ..."
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
echo "  $OUT_DIR/protein_summary.tsv"
