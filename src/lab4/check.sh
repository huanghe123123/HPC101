#!/bin/bash
# check.sh — validate Lab 4 output against golden results.
# Usage: ./check.sh [RESULT_DIR] [GOLDEN_DIR]
# See `python3 scripts/check_result.py --help` for details.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${AMSS_ROOT_DIR:-$SCRIPT_DIR}"
case "$ROOT_DIR" in
  /*) ;;
  *) ROOT_DIR="$SCRIPT_DIR/$ROOT_DIR" ;;
esac

exec "${PYTHON:-python3}" "$ROOT_DIR/scripts/check_result.py" \
    --time-tolerance "${TIME_TOLERANCE:-1e-8}" \
    "$@"
