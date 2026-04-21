#!/bin/bash
# Run MSstatsTMT pipeline → gene-level protein matrix
# Usage:
#   ./run_msstats_tmt_gene_matrix.sh [study_id] [max_runs]
# Example:
#   ./run_msstats_tmt_gene_matrix.sh PDC000120
#   ./run_msstats_tmt_gene_matrix.sh PDC000120 3   # quick run: first 3 runs only

set -e
DATA_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DATA_DIR"

STUDY="${1:-PDC000120}"
MAX_RUNS="${2:-}"
PSM_DIR="pdc_psm/${STUDY}"
OUTDIR="results/${STUDY}"

# Reference channel: 131 for TMT10, 126C for TMT11 (see sample_files_msstats_tmt.csv)
case "$STUDY" in
  PDC000120|PDC000110|PDC000116|PDC000118) REF_CHANNEL="131" ;;
  *) REF_CHANNEL="131" ;;
esac

if [[ ! -d "$PSM_DIR" ]]; then
  echo "PSM dir not found: $PSM_DIR"
  exit 1
fi

mkdir -p "$OUTDIR"

# Resolve catalog path: absolute, relative to data/, or under CPTAC_LOCAL_MIRROR/
resolve_sample_txt() {
  local raw="${1//\"/}"
  raw="${raw//$'\r'/}"
  [[ -z "$raw" ]] && return 1
  if [[ -f "$raw" ]]; then echo "$raw"; return 0; fi
  local try="$DATA_DIR/$raw"
  if [[ -f "$try" ]]; then echo "$try"; return 0; fi
  if [[ -n "${CPTAC_LOCAL_MIRROR:-}" ]]; then
    try="${CPTAC_LOCAL_MIRROR%/}/$raw"
    if [[ -f "$try" ]]; then echo "$try"; return 0; fi
  fi
  return 1
}

# Look up sample.txt path from catalog (column 2 = path, column 3 = file_name) for this study
SAMPLE_TXT=""
if [[ -f "sample_files_msstats_tmt.csv" ]]; then
  RAW_PATH=$(awk -F',' -v study="$STUDY" '$1==study {print $2; exit}' sample_files_msstats_tmt.csv)
  RAW_PATH="${RAW_PATH%"${RAW_PATH##*[![:space:]]}"}"
  if [[ -n "$RAW_PATH" ]] && resolved=$(resolve_sample_txt "$RAW_PATH"); then
    SAMPLE_TXT="$resolved"
  else
    SAMPLE_TXT=""
  fi
  # Fallback: try pdc_psm/STUDY/<file_name> or results/STUDY/<file_name>
  if [[ -z "$SAMPLE_TXT" ]]; then
    FNAME=$(awk -F',' -v study="$STUDY" '$1==study {print $3; exit}' sample_files_msstats_tmt.csv)
    FNAME="${FNAME%"${FNAME##*[![:space:]]}"}"
    for try in "${PSM_DIR}/${FNAME}" "${OUTDIR}/${FNAME}"; do
      if [[ -n "$FNAME" && -f "$try" ]]; then
        SAMPLE_TXT="$try"
        break
      fi
    done
  fi
fi
if [[ -n "$SAMPLE_TXT" ]]; then
  echo "Using study design: $SAMPLE_TXT"
fi

EXTRA=()
[[ -n "$MAX_RUNS" ]] && EXTRA+=(--max_runs "$MAX_RUNS")
[[ -n "$SAMPLE_TXT" ]] && EXTRA+=(--sample_txt "$SAMPLE_TXT")

if [[ -f "${OUTDIR}/annotation_filled.csv" ]]; then
  echo "Using existing annotation: ${OUTDIR}/annotation_filled.csv"
  Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R \
    --psm_dir "$PSM_DIR" \
    --annotation "${OUTDIR}/annotation_filled.csv" \
    --outdir "$OUTDIR" \
    --species Hs \
    "${EXTRA[@]}"
else
  echo "Auto-filling annotation (reference channel $REF_CHANNEL = Norm)"
  Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R \
    --psm_dir "$PSM_DIR" \
    --reference_channel "$REF_CHANNEL" \
    --outdir "$OUTDIR" \
    --species Hs \
    "${EXTRA[@]}"
fi

echo "Done. Outputs in $OUTDIR:"
echo "  gene_matrix.csv       - sample x gene matrix (gene-level summarization)"
echo "  protein_summary.tsv   - protein-level abundances"
echo "  qc_summary.txt        - counts and intensity distributions"
