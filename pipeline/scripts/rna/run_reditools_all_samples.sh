#!/usr/bin/env bash
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
MANIFEST="$SCRIPT_DIR/../../config/samples.tsv"
COVERAGE_HELPER="$SCRIPT_DIR/generate_reditools_coverage_limited.sh"
REBUILD_HELPER="$SCRIPT_DIR/rebuild_reditools_file_list.py"
PARALLEL="$REDITOOLS/src/cineca/parallel_reditools.py"
MERGE="$REDITOOLS/merge.sh"

SPLIT_DIR="$PROJECT/bam/splitncigarreads"
COVERAGE_ROOT="$PROJECT/reditools/coverage"
TEMP_ROOT="$PROJECT/reditools/tmp"
TABLE_DIR="$PROJECT/reditools/tables"
LOG_DIR="$PROJECT/logs"
mkdir -p "$COVERAGE_ROOT" "$TEMP_ROOT" "$TABLE_DIR" "$LOG_DIR"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

die() {
  log "ERROR: $*"
  exit 1
}

for command in samtools mpirun bgzip tabix gzip zcat python3; do
  command -v "$command" >/dev/null 2>&1 || die "Required command not found: $command"
done

[[ -f "$REFERENCE" ]] || die "Reference not found: $REFERENCE"
[[ -f "${REFERENCE}.fai" ]] || die "Reference index not found: ${REFERENCE}.fai"
[[ -x "$PYTHON_BIN" ]] || die "Python not executable: $PYTHON_BIN"
[[ -f "$PARALLEL" ]] || die "Parallel REDItools2 script not found: $PARALLEL"
[[ -f "$MERGE" ]] || die "REDItools2 merge script not found: $MERGE"
[[ -f "$REBUILD_HELPER" ]] || die "Recovery helper not found: $REBUILD_HELPER"
[[ -f "$MANIFEST" ]] || die "Sample manifest not found: $MANIFEST"

HEADER=$'Region\tPosition\tReference\tStrand\tCoverage-q30\tMeanQ\tBaseCount[A,C,G,T]\tAllSubs\tFrequency\tgCoverage-q30\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubs\tgFrequency'

output_is_valid() {
  local output="$1"
  [[ -s "$output" && -s "${output}.tbi" ]] || return 1
  gzip -t "$output" >/dev/null 2>&1 || return 1
  python3 - "$output" "$HEADER" <<'PY' >/dev/null 2>&1
import gzip
import sys
path, expected = sys.argv[1:3]
with gzip.open(path, "rt") as handle:
    observed = handle.readline().rstrip("\n")
raise SystemExit(0 if observed == expected else 1)
PY
}

coverage_is_complete() {
  local coverage_dir="$1"
  local coverage_file="$2"
  [[ -s "$coverage_file" ]] || return 1
  while IFS=$'\t' read -r chrom _rest; do
    [[ -f "${coverage_dir}${chrom}" ]] || return 1
  done < "${REFERENCE}.fai"
}

recover_file_list() {
  local temp_dir="$1"
  python3 "$REBUILD_HELPER" \
    --temp-dir "$temp_dir" \
    --fai "${REFERENCE}.fai" \
    --output "${temp_dir}files.txt"
}

merge_sample() {
  local sample="$1"
  local temp_dir="$2"
  local output="$3"
  local body="$TABLE_DIR/${sample}.body.txt.gz"
  local output_tmp="${output}.tmp.$$"

  recover_file_list "$temp_dir" \
    > "$LOG_DIR/${sample}.recover_file_list.log" 2>&1 \
    || return 1

  rm -f "$body" "${body}.tbi" "$output_tmp" "${output_tmp}.tbi"
  log "Merge REDItools2 temporary results for $sample"
  bash "$MERGE" "$temp_dir" "$body" "$THREADS" \
    > "$LOG_DIR/${sample}.merge.log" 2>&1 \
    || return 1
  [[ -s "$body" ]] || return 1

  { printf '%s\n' "$HEADER"; zcat "$body"; } | bgzip -c -@ "$THREADS" > "$output_tmp"
  gzip -t "$output_tmp" >/dev/null 2>&1 || return 1
  mv -f "$output_tmp" "$output"
  rm -f "${output}.tbi"
  tabix -f -S 1 -s 1 -b 2 -e 2 "$output"
  output_is_valid "$output"
}

while IFS=$'\t' read -r sample group replicate srr; do
  [[ "$sample" == "sample" ]] && continue

  sample_bam="$SPLIT_DIR/${sample}.splitncigarreads.bam"
  srr_bam="$SPLIT_DIR/${srr}.splitncigarreads.bam"
  if [[ -f "$sample_bam" ]]; then
    bam="$sample_bam"
  elif [[ -f "$srr_bam" ]]; then
    bam="$srr_bam"
  else
    die "Missing BAM for $sample. Checked: $sample_bam and $srr_bam"
  fi

  bam_stem="$(basename "$bam" .bam)"
  coverage_dir="$COVERAGE_ROOT/${sample}/"
  coverage_file="${coverage_dir}${bam_stem}.cov"
  coverage_marker="${coverage_dir}.complete"
  temp_dir="$TEMP_ROOT/${sample}/"
  output="$TABLE_DIR/${sample}.txt.gz"

  samtools quickcheck -v "$bam" || die "BAM integrity check failed: $bam"

  if output_is_valid "$output"; then
    log "Existing valid REDItools2 output: $sample"
    continue
  fi

  if coverage_is_complete "$coverage_dir" "$coverage_file"; then
    touch "$coverage_marker"
  else
    if [[ -d "${coverage_dir%/}" ]]; then
      archived="${coverage_dir%/}.incomplete.$(date '+%Y%m%d_%H%M%S')"
      log "Archive incomplete coverage directory: $archived"
      mv "${coverage_dir%/}" "$archived"
    fi
    log "Coverage: $sample"
    bash "$COVERAGE_HELPER" "$bam" "$coverage_dir" "${REFERENCE}.fai" "$COVERAGE_JOBS" \
      > "$LOG_DIR/${sample}.coverage.log" 2>&1
    coverage_is_complete "$coverage_dir" "$coverage_file" \
      || die "Coverage generation was incomplete for $sample"
    touch "$coverage_marker"
  fi

  if [[ -d "$temp_dir" ]] && merge_sample "$sample" "$temp_dir" "$output"; then
    log "Recovered completed REDItools2 temporary output: $sample"
    continue
  fi

  if [[ -d "${temp_dir%/}" ]]; then
    archived="${temp_dir%/}.incomplete.$(date '+%Y%m%d_%H%M%S')"
    log "Archive incomplete REDItools2 temporary directory: $archived"
    mv "${temp_dir%/}" "$archived"
  fi
  mkdir -p "$temp_dir"
  rm -f "$output" "${output}.tbi"

  log "REDItools2: $sample"
  set +e
  mpirun -np "$MPI_PROCS" "$PYTHON_BIN" "$PARALLEL" \
    -f "$bam" -r "$REFERENCE" \
    -G "$coverage_file" -D "$coverage_dir" \
    -t "$temp_dir" -Z "${REFERENCE}.fai" \
    -S -me 20 \
    > "$LOG_DIR/${sample}.REDItools2.log" 2>&1
  mpi_status=$?
  set -e

  if merge_sample "$sample" "$temp_dir" "$output"; then
    if [[ $mpi_status -ne 0 ]]; then
      log "Recovered $sample after REDItools2 exited with status $mpi_status"
    fi
    log "Completed: $output"
    continue
  fi

  die "REDItools2 did not produce a complete recoverable result for $sample (MPI status $mpi_status). See $LOG_DIR/${sample}.REDItools2.log"
done < "$MANIFEST"
