set -Eeuo pipefail

if [[ $# -lt 5 || $# -gt 7 ]]; then
  echo "Usage: bash $0 PROJECT REFERENCE REDITOOLS PYTHON MPI_PROCS [COVERAGE_JOBS] [THREADS]" >&2
  exit 2
fi

PROJECT="$(readlink -m "$1")"
REFERENCE="$(readlink -m "$2")"
REDITOOLS="$(readlink -m "$3")"
PYTHON_BIN="$(readlink -m "$4")"
MPI_PROCS="$5"
COVERAGE_JOBS="${6:-8}"
THREADS="${7:-8}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="$SCRIPT_DIR/../config/samples.tsv"
COVERAGE_HELPER="$SCRIPT_DIR/generate_reditools_coverage_limited.sh"
PARALLEL="$REDITOOLS/src/cineca/parallel_reditools.py"
MERGE="$REDITOOLS/merge.sh"

SPLIT_DIR="$PROJECT/bam/splitncigarreads"
COVERAGE_ROOT="$PROJECT/reditools/coverage"
TEMP_ROOT="$PROJECT/reditools/tmp"
TABLE_DIR="$PROJECT/reditools/tables"
LOG_DIR="$PROJECT/logs"
mkdir -p "$COVERAGE_ROOT" "$TEMP_ROOT" "$TABLE_DIR" "$LOG_DIR"

[[ -f "$REFERENCE" ]] || { echo "Reference not found: $REFERENCE" >&2; exit 1; }
[[ -f "${REFERENCE}.fai" ]] || { echo "Reference index not found" >&2; exit 1; }
[[ -x "$PYTHON_BIN" ]] || { echo "Python not executable: $PYTHON_BIN" >&2; exit 1; }
[[ -f "$PARALLEL" ]] || { echo "Parallel REDItools2 script not found" >&2; exit 1; }
[[ -f "$MERGE" ]] || { echo "REDItools2 merge script not found" >&2; exit 1; }

HEADER=$'Region\tPosition\tReference\tStrand\tCoverage-q30\tMeanQ\tBaseCount[A,C,G,T]\tAllSubs\tFrequency\tgCoverage-q30\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubs\tgFrequency'

while IFS=$'\t' read -r sample group replicate srr; do
  [[ "$sample" == "sample" ]] && continue

  sample_bam="$SPLIT_DIR/${sample}.splitncigarreads.bam"
  srr_bam="$SPLIT_DIR/${srr}.splitncigarreads.bam"
  if [[ -f "$sample_bam" ]]; then
    bam="$sample_bam"
  elif [[ -f "$srr_bam" ]]; then
    bam="$srr_bam"
  else
    echo "Missing BAM for $sample. Checked: $sample_bam and $srr_bam" >&2
    exit 1
  fi

  bam_stem="$(basename "$bam" .bam)"
  coverage_dir="$COVERAGE_ROOT/${sample}/"
  coverage_file="${coverage_dir}${bam_stem}.cov"
  temp_dir="$TEMP_ROOT/${sample}"
  output="$TABLE_DIR/${sample}.txt.gz"

  samtools quickcheck -v "$bam"

  if [[ ! -s "$coverage_file" ]]; then
    echo "[$(date '+%F %T')] Coverage: $sample" >&2
    bash "$COVERAGE_HELPER" "$bam" "$coverage_dir" "${REFERENCE}.fai" "$COVERAGE_JOBS" \
      > "$LOG_DIR/${sample}.coverage.log" 2>&1
  fi
  [[ -s "$coverage_file" ]] || { echo "Missing coverage file: $coverage_file" >&2; exit 1; }

  if [[ -s "$output" && -s "${output}.tbi" ]]; then
    echo "[$(date '+%F %T')] Existing REDItools2 output: $sample" >&2
    continue
  fi

  mkdir -p "$temp_dir"
  echo "[$(date '+%F %T')] REDItools2: $sample" >&2
  mpirun -np "$MPI_PROCS" "$PYTHON_BIN" "$PARALLEL" \
    -f "$bam" -r "$REFERENCE" \
    -G "$coverage_file" -D "$coverage_dir" \
    -t "$temp_dir" -Z "${REFERENCE}.fai" \
    -S -me 20 \
    > "$LOG_DIR/${sample}.REDItools2.log" 2>&1

  body="$TABLE_DIR/${sample}.body.txt.gz"
  bash "$MERGE" "$temp_dir" "$body" "$THREADS" \
    > "$LOG_DIR/${sample}.merge.log" 2>&1
  [[ -s "$body" ]] || { echo "Empty merged table: $body" >&2; exit 1; }

  { printf '%s\n' "$HEADER"; zcat "$body"; } | bgzip -c -@ "$THREADS" > "$output"
  tabix -f -S 1 -s 1 -b 2 -e 2 "$output"
  echo "[$(date '+%F %T')] Completed: $output" >&2
done < "$MANIFEST"
