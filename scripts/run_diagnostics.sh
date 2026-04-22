#!/usr/bin/env bash
# Subset of benchmark diagnostics (fast checks). Full suite is inside run_overnight_v2.sh.
#
# Usage:
#   ./scripts/run_diagnostics.sh preflight   # gene coverage + intersection lists (R)
#   ./scripts/run_diagnostics.sh structure   # Python structure batch (PCA-related summaries)
#   ./scripts/run_diagnostics.sh all         # preflight then structure
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
export PYTHONPATH="${REPO}/src${PYTHONPATH:+:$PYTHONPATH}"
if [[ -x "${REPO}/.venv/bin/python3" ]]; then
  PYTHON="${REPO}/.venv/bin/python3"
else
  PYTHON=python3
fi

UNION_DIR="${PROCESSED_UNION_DIR:-${REPO}/data/processed/union}"
INTER_DIR="${INTERSECTION_DIR:-${REPO}/data/processed}"

preflight() {
  Rscript "${REPO}/scripts/benchmark/preflight_diagnostics.R" --repo-root "$REPO" \
    --processed-dir "$UNION_DIR" || return 1
  Rscript "${REPO}/scripts/benchmark/compute_intersection_masks.R" --repo-root "$REPO" \
    --intersection-out-dir "$UNION_DIR" || return 1
}

structure() {
  "$PYTHON" "${REPO}/scripts/benchmark/run_structure_batch.py" --repo "$REPO" \
    --meta-dir "${META_DIR:-$UNION_DIR}" || return 1
}

case "${1:-all}" in
  preflight) preflight ;;
  structure) structure ;;
  all) preflight && structure ;;
  *)
    echo "Usage: $0 [preflight|structure|all]"
    exit 1
    ;;
esac
