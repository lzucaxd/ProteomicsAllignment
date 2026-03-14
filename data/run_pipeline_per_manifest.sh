#!/usr/bin/env bash
# Run pipeline one manifest at a time: download PSM, then MSstatsTMT → gene matrix.
# Uses manifests in ./manifests/ (no duplicates). Output: pdc_psm/{study_id}/, results/{study_id}/
# Requires: Python with 'requests' (pip install requests), R with MSstatsTMT.
#
# Optional (disk-efficient run):
#   --cleanup-after   After each study succeeds, remove heavy intermediates (parsed_psm_long.tsv, msstats_input.tsv).
#   --delete-psm     When used with --cleanup-after, also remove pdc_psm/{study_id}/ after the pipeline.
# Example: ./run_pipeline_per_manifest.sh --cleanup-after --delete-psm

set -e
DATA_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DATA_DIR"
MANIFESTS_DIR="$DATA_DIR/manifests"
PSM_ROOT="$DATA_DIR/pdc_psm"
RESULTS_ROOT="$DATA_DIR/results"

CLEANUP_AFTER=false
DELETE_PSM=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cleanup-after) CLEANUP_AFTER=true; shift ;;
    --delete-psm)    DELETE_PSM=true; shift ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Prefer venv Python if it has requests
if [ -x "$DATA_DIR/.venv/bin/python" ] && "$DATA_DIR/.venv/bin/python" -c "import requests" 2>/dev/null; then
  PYTHON="$DATA_DIR/.venv/bin/python"
else
  PYTHON=python3
fi

mkdir -p "$PSM_ROOT" "$RESULTS_ROOT"

# Get study ID from first data row of a manifest CSV
get_study_id() {
  local manifest="$1"
  awk -F',' 'NR==1 { for(i=1;i<=NF;i++) if($i~/PDC Study ID/) col=i; next } col && $col!="" { gsub(/"/,"",$col); print $col; exit }' "$manifest"
}

for manifest in "$MANIFESTS_DIR"/PDC_file_manifest_*.csv; do
  [ -f "$manifest" ] || continue
  name="$(basename "$manifest")"
  study_id="$(get_study_id "$manifest")"
  if [ -z "$study_id" ]; then
    echo "Skip $name: could not get PDC Study ID"
    continue
  fi
  echo "=============================================="
  echo "Manifest: $name  →  Study: $study_id"
  echo "=============================================="

  echo "[1/2] Downloading PSM files..."
  "$PYTHON" pdc_manifest_downloader.py \
    --manifest "$manifest" \
    --outdir "$PSM_ROOT" \
    --include-category "Peptide Spectral Matches" \
    --ext .psm

  psm_dir="$PSM_ROOT/$study_id"
  outdir="$RESULTS_ROOT/$study_id"
  if [ ! -d "$psm_dir" ]; then
    echo "No PSM dir after download: $psm_dir (no .psm in this manifest?)"
    continue
  fi

  echo "[check] Sample file for MSstatsTMT..."
  if ! "$PYTHON" "$DATA_DIR/check_studies_sample_file.py" --study "$study_id" 2>/dev/null; then
    echo "  WARNING: $study_id has NO entry in sample_files_msstats_tmt.csv — pipeline may need --reference_channel or manual annotation."
  fi

  echo "[2/2] Running MSstatsTMT pipeline..."
  Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R \
    --psm_dir "$psm_dir" \
    --outdir "$outdir"

  echo "Done: $study_id"

  if [ "$CLEANUP_AFTER" = true ] && [ -f "$outdir/gene_matrix.csv" ]; then
    echo "[cleanup] Freeing disk for $study_id..."
    rm -f "$outdir/parsed_psm_long.tsv" "$outdir/msstats_input.tsv"
    $DELETE_PSM && [ -d "$psm_dir" ] && rm -rf "$psm_dir" && echo "  removed pdc_psm/$study_id"
    echo "  done."
  fi
  echo ""
done

echo "All manifests processed."
