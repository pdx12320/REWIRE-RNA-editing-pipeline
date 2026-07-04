set -Eeuo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: bash $0 BAM COVERAGE_DIR REFERENCE.fai [JOBS]" >&2
  exit 2
fi

BAM="$(readlink -m "$1")"
COVERAGE_DIR="$(readlink -m "$2")/"
FAI="$(readlink -m "$3")"
JOBS="${4:-8}"
SAMTOOLS_BIN="$(command -v samtools)"

[[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
[[ -f "$FAI" ]] || { echo "ERROR: FAI not found: $FAI" >&2; exit 1; }
[[ -n "$SAMTOOLS_BIN" ]] || { echo "ERROR: samtools not found" >&2; exit 1; }
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: JOBS must be a positive integer" >&2; exit 1; }

samtools quickcheck -v "$BAM" || { echo "ERROR: BAM integrity check failed: $BAM" >&2; exit 1; }

FILE_ID="$(basename "$BAM")"
FILE_ID="${FILE_ID%.bam}"
FINAL_COV="${COVERAGE_DIR}${FILE_ID}.cov"
TEMP_COV="${FINAL_COV}.tmp.$$"
COMPLETE_MARKER="${COVERAGE_DIR}.complete"

mkdir -p "$COVERAGE_DIR"
rm -f "$COMPLETE_MARKER" "$TEMP_COV"
export BAM COVERAGE_DIR SAMTOOLS_BIN

cut -f1 "$FAI" | xargs -r -n 1 -P "$JOBS" bash -c '
  chrom="$1"
  printf "[%s] Coverage %s\n" "$(date "+%F %T")" "$chrom" >&2
  temporary="${COVERAGE_DIR}.${chrom}.tmp.$$"
  "$SAMTOOLS_BIN" depth -r "$chrom" "$BAM" > "$temporary"
  mv -f "$temporary" "${COVERAGE_DIR}${chrom}"
' _

: > "$TEMP_COV"
while IFS=$'\t' read -r chrom _rest; do
  part="${COVERAGE_DIR}${chrom}"
  [[ -f "$part" ]] || { echo "ERROR: missing per-contig file: $part" >&2; exit 1; }
  cat "$part" >> "$TEMP_COV"
done < "$FAI"

[[ -s "$TEMP_COV" ]] || { echo "ERROR: final coverage file is empty: $TEMP_COV" >&2; exit 1; }
mv -f "$TEMP_COV" "$FINAL_COV"
touch "$COMPLETE_MARKER"
echo "Coverage complete: $FINAL_COV" >&2
