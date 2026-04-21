#!/usr/bin/env bash
# Build union matrices under data/processed/union/ after v2 annotation processing.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
if [[ ! -f "${REPO}/renv/activate.R" ]]; then
  export R_PROFILE_USER="/dev/null"
fi
export PYTHONPATH="${REPO}/src${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "${REPO}/data/processed/union" "${REPO}/reports/benchmark_master/diagnostics"

echo "[v2] Process CCLE subtype annotations"
Rscript "${REPO}/scripts/benchmark/process_ccle_annotations_v2.R" --repo-root "$REPO"

echo "[v2] Run preprocessing -> data/processed/union"
python3 "${REPO}/scripts/run_preprocessing.py" --repo-root "$REPO" --output-dir data/processed/union

echo "[v2] Preflight + intersection gene lists"
Rscript "${REPO}/scripts/benchmark/preflight_diagnostics.R" --repo-root "$REPO" \
  --processed-dir data/processed/union
Rscript "${REPO}/scripts/benchmark/compute_intersection_masks.R" --repo-root "$REPO" \
  --intersection-out-dir data/processed/union

echo "Done. Union matrices: ${REPO}/data/processed/union/"
echo "Intersection lists: ${REPO}/data/processed/intersection_genes_*.txt"
