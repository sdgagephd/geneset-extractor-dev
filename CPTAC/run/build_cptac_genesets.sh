#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
exec "$PYTHON_BIN" "$HERE/../src/run_cptac_ptm_model.py" "$@"
