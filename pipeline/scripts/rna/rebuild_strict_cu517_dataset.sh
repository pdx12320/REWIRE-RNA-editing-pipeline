#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -lt 5 || $# -gt 6 ]]; then
  echo "Usage: bash $0 PROJECT GRCH38_FASTA GENCODE_GTF EGFP_GC_REPORTER_FASTA MANIFEST [THREADS]" >&2
  exit 2
fi

PROJECT="$(readlink -m "$1")"
GENOME_FASTA="$(readlink -m "$2")"
GENOME_GTF="$(readlink -m "$3")"
REPORTER_FASTA="$(readlink -m "$4")"
MANIFEST="$(readlink -m "$5")"
THREADS="${6:-32}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for path in "$GENOME_FASTA" "$GENOME_GTF" "$REPORTER_FASTA" "$MANIFEST"; do
  [[ -f "$path" ]] || { echo "ERROR: required input not found: $path" >&2; exit 1; }
done

REFERENCE_DIR="$PROJECT/reference/strict_augmented"
STAR_INDEX="$REFERENCE_DIR/star_index"
AUGMENTED_FASTA="$REFERENCE_DIR/GRCh38_plus_EGFP_GC.fa"
AUGMENTED_GTF="$REFERENCE_DIR/GRCh38_plus_EGFP_GC.gtf"
STRICT_DIR="$PROJECT/strict_lamar"
FULL_TABLE="$STRICT_DIR/full_coverage_cytidines.tsv.gz"
MAPPABILITY_TABLE="$STRICT_DIR/mappability_audit.tsv.gz"

mkdir -p "$REFERENCE_DIR" "$STAR_INDEX" "$STRICT_DIR"

python3 "$SCRIPT_DIR/build_augmented_reference.py" \
  --genome-fasta "$GENOME_FASTA" \
  --genome-gtf "$GENOME_GTF" \
  --reporter-fasta "$REPORTER_FASTA" \
  --output-fasta "$AUGMENTED_FASTA" \
  --output-gtf "$AUGMENTED_GTF" \
  --contig-name EGFP_GC_reporter \
  --target-position 459 \
  --expected-context GC

samtools faidx "$AUGMENTED_FASTA"
DICT="${AUGMENTED_FASTA%.*}.dict"
if [[ ! -s "$DICT" ]]; then
  gatk CreateSequenceDictionary -R "$AUGMENTED_FASTA" -O "$DICT"
fi

STAR \
  --runMode genomeGenerate \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --genomeFastaFiles "$AUGMENTED_FASTA" \
  --sjdbGTFfile "$AUGMENTED_GTF" \
  --sjdbOverhang 149

python3 "$SCRIPT_DIR/run_star_alignment.py" \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest "$MANIFEST" \
  --threads "$THREADS" \
  --expected-contig EGFP_GC_reporter \
  --expected-samples 6 \
  --force

python3 "$SCRIPT_DIR/run_gatk_preprocessing.py" \
  --project "$PROJECT" \
  --reference "$AUGMENTED_FASTA" \
  --manifest "$MANIFEST" \
  --expected-contig EGFP_GC_reporter \
  --expected-samples 6 \
  --force

python3 "$SCRIPT_DIR/build_full_coverage_cytidine_labels.py" \
  --manifest "$MANIFEST" \
  --bam-dir "$PROJECT/bam/splitncigarreads" \
  --reference "$AUGMENTED_FASTA" \
  --gtf "$AUGMENTED_GTF" \
  --output "$FULL_TABLE" \
  --min-depth 50 \
  --min-mapq 30 \
  --min-baseq 20 \
  --positive-threshold 0.10 \
  --min-positive-treated-reps 2 \
  --min-positive-alt-count 3 \
  --max-positive-control-median 0.01

if [[ ! -s "${AUGMENTED_FASTA}.bwt" ]]; then
  bwa index "$AUGMENTED_FASTA"
fi

python3 "$SCRIPT_DIR/audit_sequence_mappability.py" \
  --input "$FULL_TABLE" \
  --reference "$AUGMENTED_FASTA" \
  --output "$MAPPABILITY_TABLE" \
  --sequence-column sequence_context \
  --eligible-only \
  --threads "$THREADS"

python3 "$SCRIPT_DIR/finalize_strict_cytidine_dataset.py" \
  --mappability-table "$MAPPABILITY_TABLE" \
  --output-dir "$STRICT_DIR/dataset_1to200" \
  --negative-ratio 200 \
  --seed 20260822 \
  --overwrite

echo "Strict CU5.17 dataset complete: $STRICT_DIR/dataset_1to200" >&2
