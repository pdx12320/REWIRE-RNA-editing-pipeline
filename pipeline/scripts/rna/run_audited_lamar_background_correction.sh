#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: bash $0 PROJECT REFERENCE_FASTA [SAMTOOLS] [PYTHON] [EXTRA_ARGS...]" >&2
  exit 2
fi

PROJECT="$(readlink -m "$1")"
REFERENCE="$(readlink -m "$2")"
shift 2

SAMTOOLS="${1:-$(command -v samtools || true)}"
[[ $# -gt 0 ]] && shift
PYTHON="${1:-$(command -v python3 || true)}"
[[ $# -gt 0 ]] && shift

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
SCRIPT="$SCRIPT_DIR/run_audited_lamar_background_correction.py"
TEST_FILE="$REPO_ROOT/tests/test_audited_lamar_background_correction.py"

[[ -d "$PROJECT" ]] || { echo "ERROR: project directory not found: $PROJECT" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "ERROR: reference FASTA not found: $REFERENCE" >&2; exit 1; }
[[ -x "$SAMTOOLS" ]] || { echo "ERROR: samtools executable not found: $SAMTOOLS" >&2; exit 1; }
[[ -x "$PYTHON" ]] || { echo "ERROR: Python executable not found: $PYTHON" >&2; exit 1; }
[[ -f "$TEST_FILE" ]] || { echo "ERROR: audit tests not found: $TEST_FILE" >&2; exit 1; }

exec "$PYTHON" "$SCRIPT" \
  --project "$PROJECT" \
  --reference "$REFERENCE" \
  --samtools "$SAMTOOLS" \
  --python "$PYTHON" \
  --test-file "$TEST_FILE" \
  --min-mapq 30 \
  --min-baseq 20 \
  --min-coverage 20 \
  --min-group-replicates 2 \
  --sequence-length 101 \
  --threads 4 \
  "$@"
