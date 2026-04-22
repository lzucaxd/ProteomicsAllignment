#!/usr/bin/env bash
# Harmonization method smoke test / registry check (Python harmonize package).
# For overnight **union** method matrices (raw, bridge_shift, bridge_scale, celligner), use:
#   python scripts/benchmark/regenerate_methods_union.py
# after union preprocessing outputs exist (see scripts/benchmark/run_overnight_v2.sh Step 2).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
export PYTHONPATH="${REPO}/src${PYTHONPATH:+:$PYTHONPATH}"
if [[ -x "${REPO}/.venv/bin/python3" ]]; then
  exec "${REPO}/.venv/bin/python3" "${REPO}/scripts/run_methods.py" "$@"
fi
exec python3 "${REPO}/scripts/run_methods.py" "$@"
