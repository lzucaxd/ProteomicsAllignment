#!/usr/bin/env bash
# Optional: free disk space after pipeline has produced gene_matrix.csv for a study.
# Removes heavy intermediates (parsed_psm_long.tsv, msstats_input.tsv) and optionally
# the raw PSM download for the study. Does NOT remove gene_matrix.csv or annotation.
#
# Usage:
#   ./cleanup_study_disk.sh PDC000120              # one study: intermediates only
#   ./cleanup_study_disk.sh PDC000120 --delete-psm  # one study: intermediates + pdc_psm/PDC000120
#   ./cleanup_study_disk.sh --all                   # all studies in results/: intermediates only
#   ./cleanup_study_disk.sh --all --delete-psm      # all studies: intermediates + pdc_psm
#
# Requires: gene_matrix.csv present in results/{study_id}/ (otherwise that study is skipped).

set -e
DATA_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DATA_DIR"
RESULTS_ROOT="$DATA_DIR/results"
PSM_ROOT="$DATA_DIR/pdc_psm"

DELETE_PSM=false
STUDIES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all)
      for d in "$RESULTS_ROOT"/PDC*/; do
        [ -d "$d" ] || continue
        id=$(basename "$d")
        [ -f "$RESULTS_ROOT/$id/gene_matrix.csv" ] && STUDIES+=("$id")
      done
      shift
      ;;
    --delete-psm)
      DELETE_PSM=true
      shift
      ;;
    PDC*)
      STUDIES+=("$1")
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

if [[ ${#STUDIES[@]} -eq 0 ]]; then
  echo "Usage: $0 [--all] [--delete-psm] [PDC000120 ...]" >&2
  echo "  --all         cleanup all studies in results/ that have gene_matrix.csv" >&2
  echo "  --delete-psm  also remove pdc_psm/{study_id}/ for each study" >&2
  exit 1
fi

for study in "${STUDIES[@]}"; do
  outdir="$RESULTS_ROOT/$study"
  if [ ! -f "$outdir/gene_matrix.csv" ]; then
    echo "Skip $study: no gene_matrix.csv"
    continue
  fi
  echo "Cleanup $study ..."
  [ -f "$outdir/parsed_psm_long.tsv" ] && rm -f "$outdir/parsed_psm_long.tsv" && echo "  removed parsed_psm_long.tsv"
  [ -f "$outdir/msstats_input.tsv" ]   && rm -f "$outdir/msstats_input.tsv"   && echo "  removed msstats_input.tsv"
  if $DELETE_PSM && [ -d "$PSM_ROOT/$study" ]; then
    rm -rf "$PSM_ROOT/$study"
    echo "  removed pdc_psm/$study"
  fi
  echo "  done."
done
