#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SQANTI3_QC_SCRIPT="$1"
shift
source "$SCRIPT_DIR/sqanti3_env.sh"
exec python "$SCRIPT_DIR/run_sqanti3.py" "$SQANTI3_QC_SCRIPT" "$@"
