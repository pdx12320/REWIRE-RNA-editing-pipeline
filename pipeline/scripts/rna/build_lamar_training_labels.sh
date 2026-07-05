#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: bash $0 PROJECT REFERENCE_FASTA CANDIDATE_TABLE MANIFEST [SITE_METADATA]" >&2
  exit 2
fi

PROJECT="$(readlink -m "$1")"
REFERENCE="$(readlink -m "$2")"
CANDIDATE_TABLE="$(readlink -m "$3")"
MANIFEST="$(readlink -m "$4")"
SITE_METADATA="${5:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAM_DIR="$PROJECT/bam/splitncigarreads"
COUNT_DIR="$PROJECT/lamar_training/base_counts"
OUTPUT="$PROJECT/lamar_training/CU5.17_EGFP_GC.lamar_training_labels.tsv.gz"

for path in "$REFERENCE" "$CANDIDATE_TABLE" "$MANIFEST"; do
  [[ -f "$path" ]] || { echo "ERROR: required input not found: $path" >&2; exit 1; }
done
[[ -d "$BAM_DIR" ]] || { echo "ERROR: BAM directory not found: $BAM_DIR" >&2; exit 1; }
mkdir -p "$COUNT_DIR" "$(dirname "$OUTPUT")"

python3 "$SCRIPT_DIR/pileup_candidate_base_counts.py" \
  --manifest "$MANIFEST" \
  --candidate-table "$CANDIDATE_TABLE" \
  --bam-dir "$BAM_DIR" \
  --reference "$REFERENCE" \
  --output-dir "$COUNT_DIR" \
  --min-mapq 30 \
  --min-baseq 20

build_args=(
  --manifest "$MANIFEST"
  --candidate-table "$CANDIDATE_TABLE"
  --count-dir "$COUNT_DIR"
  --reference "$REFERENCE"
  --output "$OUTPUT"
  --window-size 101
  --min-allele-depth 20
)
if [[ -n "$SITE_METADATA" ]]; then
  SITE_METADATA="$(readlink -m "$SITE_METADATA")"
  [[ -f "$SITE_METADATA" ]] || { echo "ERROR: site metadata not found: $SITE_METADATA" >&2; exit 1; }
  build_args+=(--site-metadata "$SITE_METADATA")
fi

python3 "$SCRIPT_DIR/build_lamar_training_table.py" "${build_args[@]}"
echo "LAMAR training labels complete: $OUTPUT" >&2
