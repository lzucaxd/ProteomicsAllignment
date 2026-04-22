#!/usr/bin/env bash
# Canonical full benchmark (union matrices, methods, limma, nulls, ceilings, tables).
# Wrapper only — implementation lives in scripts/benchmark/run_overnight_v2.sh
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec bash "${REPO}/scripts/benchmark/run_overnight_v2.sh" "$@"
