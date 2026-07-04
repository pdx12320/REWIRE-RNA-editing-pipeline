#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
cat <<'USAGE'
Usage:
  bash scripts/wgs/run_3run_wgs_pipeline.sh \
    --reference /path/GRCh38.fa \
    --outdir /path/hek293t_wgs \
    [--runs config/wgs_runs.tsv] \
    [--mode auto|merge|consensus] \
    [--threads 32]

Modes:
  auto       Query ENA metadata. Merge only when all runs share one BioSample;
             otherwise call each run and build a >=2/3 consensus blacklist.
  merge      Treat all runs as sequencing runs from one biological sample.
  consensus  Treat runs as independent public 293T genomes and retain exact
             SNVs supported by at least 2 runs.
USAGE
}

RUNS_FILE="config/wgs_runs.tsv"
REFERENCE=""
OUTDIR=""
MODE="auto"
THREADS=32
JAVA_OPTIONS="-Xmx24g"
MIN_DP=10
MIN_ALT=3
MIN_VAF=0.05
MIN_QUAL=20

while [[ $# -gt 0 ]]; do
  case "$1" in
    --runs) RUNS_FILE="$2"; shift 2;;
    --reference) REFERENCE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --java-options) JAVA_OPTIONS="$2"; shift 2;;
    --min-dp) MIN_DP="$2"; shift 2;;
    --min-alt) MIN_ALT="$2"; shift 2;;
    --min-vaf) MIN_VAF="$2"; shift 2;;
    --min-qual) MIN_QUAL="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "ERROR: unknown argument $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$REFERENCE" && -f "$REFERENCE" ]] || { echo "ERROR: valid --reference is required" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 1; }
[[ "$MODE" =~ ^(auto|merge|consensus)$ ]] || { echo "ERROR: invalid --mode" >&2; exit 1; }

for cmd in prefetch fasterq-dump vdb-validate pigz samtools bcftools gatk python3 curl bgzip tabix; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing command: $cmd" >&2; exit 1; }
done
if command -v bwa-mem2 >/dev/null 2>&1; then
  BWA=(bwa-mem2)
elif command -v bwa >/dev/null 2>&1; then
  BWA=(bwa)
else
  echo "ERROR: bwa-mem2 or bwa is required" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_FILE="$(readlink -m "$RUNS_FILE")"
REFERENCE="$(readlink -m "$REFERENCE")"
OUTDIR="$(readlink -m "$OUTDIR")"
mkdir -p "$OUTDIR"/{metadata,sra,fastq,bam,vcf,qc,logs,tmp}

[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"
DICT="${REFERENCE%.*}.dict"
[[ -f "$DICT" ]] || gatk --java-options "$JAVA_OPTIONS" CreateSequenceDictionary -R "$REFERENCE" -O "$DICT"

META="$OUTDIR/metadata/sra_metadata.tsv"
bash "$SCRIPT_DIR/00_check_sra_metadata.sh" "$RUNS_FILE" "$META" | tee "$OUTDIR/logs/metadata.log"

if [[ "$MODE" == "auto" ]]; then
  MODE="$(python3 - "$META" <<'PY'
import csv,sys
with open(sys.argv[1]) as fh:
    rows=list(csv.DictReader(fh,delimiter='\t'))
samples={r.get('sample_accession','') for r in rows if r.get('sample_accession','')}
strategies={r.get('library_strategy','').upper() for r in rows}
layouts={r.get('library_layout','').upper() for r in rows}
if strategies != {'WGS'} or layouts != {'PAIRED'}:
    raise SystemExit('Metadata check failed: all runs must be paired-end WGS')
print('merge' if len(samples)==1 and len(rows)==3 else 'consensus')
PY
)"
fi
echo "Resolved mode: $MODE" | tee "$OUTDIR/logs/mode.log"

mapfile -t RUNS < <(awk -F '\t' 'NR>1 && $1!="" {print $1}' "$RUNS_FILE")
[[ ${#RUNS[@]} -gt 0 ]] || { echo "ERROR: no runs in $RUNS_FILE" >&2; exit 1; }

declare -A SAMPLE_ACC EXP_ACC
while IFS=$'\t' read -r run study sample experiment rest; do
  [[ "$run" == "run_accession" ]] && continue
  SAMPLE_ACC["$run"]="$sample"
  EXP_ACC["$run"]="$experiment"
done < "$META"

if [[ "${BWA[0]}" == "bwa-mem2" ]]; then
  [[ -f "${REFERENCE}.0123" ]] || bwa-mem2 index "$REFERENCE"
else
  [[ -f "${REFERENCE}.bwt" ]] || bwa index "$REFERENCE"
fi

ALIGNED_BAMS=()
for run in "${RUNS[@]}"; do
  echo "[$(date '+%F %T')] Processing $run"
  RUN_DIR="$OUTDIR/sra/$run"
  mkdir -p "$RUN_DIR" "$OUTDIR/tmp/$run"

  if [[ ! -e "$RUN_DIR/$run.sra" ]]; then
    prefetch "$run" --max-size u -O "$RUN_DIR" 2>&1 | tee "$OUTDIR/logs/${run}.prefetch.log"
  fi
  SRA_PATH="$(find "$RUN_DIR" -type f -name '*.sra' -print -quit)"
  [[ -n "$SRA_PATH" ]] || { echo "ERROR: .sra file not found for $run" >&2; exit 1; }
  vdb-validate "$SRA_PATH" 2>&1 | tee "$OUTDIR/logs/${run}.validate.log"

  R1="$OUTDIR/fastq/${run}_1.fastq.gz"
  R2="$OUTDIR/fastq/${run}_2.fastq.gz"
  if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    fasterq-dump "$SRA_PATH" --split-files --threads "$THREADS" \
      --outdir "$OUTDIR/fastq" --temp "$OUTDIR/tmp/$run" \
      2>&1 | tee "$OUTDIR/logs/${run}.fasterq.log"
    [[ -s "$OUTDIR/fastq/${run}_1.fastq" && -s "$OUTDIR/fastq/${run}_2.fastq" ]] || {
      echo "ERROR: paired FASTQ files were not produced for $run" >&2
      exit 1
    }
    pigz -p "$THREADS" "$OUTDIR/fastq/${run}_1.fastq" "$OUTDIR/fastq/${run}_2.fastq"
  fi

  RAW_BAM="$OUTDIR/bam/${run}.sorted.bam"
  if [[ ! -s "$RAW_BAM" ]]; then
    if [[ "$MODE" == "merge" ]]; then
      SM="HEK293T_PUBLIC_WGS"
    else
      SM="$run"
    fi
    LB="${EXP_ACC[$run]:-$run}"
    RG="@RG\tID:${run}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA\tPU:${run}"
    "${BWA[@]}" mem -t "$THREADS" -R "$RG" "$REFERENCE" "$R1" "$R2" \
      2> "$OUTDIR/logs/${run}.bwa.log" | \
      samtools sort -@ "$THREADS" -m 2G -o "$RAW_BAM" -
    samtools index -@ "$THREADS" "$RAW_BAM"
    samtools quickcheck -v "$RAW_BAM"
  fi
  ALIGNED_BAMS+=("$RAW_BAM")
done

call_one() {
  local name="$1" bam="$2"
  local markdup="$OUTDIR/bam/${name}.markdup.bam"
  local metrics="$OUTDIR/qc/${name}.markduplicates.metrics.txt"
  local norm="$OUTDIR/vcf/${name}.raw.norm.vcf.gz"
  local plain="$OUTDIR/vcf/${name}.filtered.SNV.vcf"
  local final="$OUTDIR/vcf/${name}.filtered.SNV.vcf.gz"

  if [[ ! -s "$markdup" ]]; then
    gatk --java-options "$JAVA_OPTIONS" MarkDuplicates \
      -I "$bam" -O "$markdup" -M "$metrics" \
      --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
      --TMP_DIR "$OUTDIR/tmp"
  fi
  samtools quickcheck -v "$markdup"

  if [[ ! -s "$norm" ]]; then
    bcftools mpileup -Ou -f "$REFERENCE" -q 20 -Q 20 -d 100000 \
      -a FORMAT/AD,FORMAT/DP "$markdup" | \
    bcftools call -mv -Ou | \
    bcftools norm -f "$REFERENCE" -m -any -Oz -o "$norm"
    bcftools index -f -t "$norm"
  fi

  python3 "$SCRIPT_DIR/filter_single_sample_vcf.py" \
    --input "$norm" --output "$plain" \
    --min-dp "$MIN_DP" --min-alt "$MIN_ALT" \
    --min-vaf "$MIN_VAF" --min-qual "$MIN_QUAL"
  bgzip -f -@ "$THREADS" "$plain"
  mv "${plain}.gz" "$final"
  tabix -f -p vcf "$final"
  echo "$final"
}

if [[ "$MODE" == "merge" ]]; then
  MERGED="$OUTDIR/bam/HEK293T_3runs.merged.bam"
  if [[ ! -s "$MERGED" ]]; then
    samtools merge -@ "$THREADS" -f "$MERGED" "${ALIGNED_BAMS[@]}"
    samtools index -@ "$THREADS" "$MERGED"
  fi
  FINAL_VCF="$(call_one HEK293T_3runs "$MERGED" | tail -n 1)"
  echo "Ready-to-use public WGS blacklist: $FINAL_VCF"
else
  FILTERED=()
  for i in "${!RUNS[@]}"; do
    FILTERED+=("$(call_one "${RUNS[$i]}" "${ALIGNED_BAMS[$i]}" | tail -n 1)")
  done
  UNION_PLAIN="$OUTDIR/vcf/HEK293T_3runs.union.SNV.vcf"
  CONS_PLAIN="$OUTDIR/vcf/HEK293T_3runs.consensus2of3.SNV.vcf"
  args=()
  for v in "${FILTERED[@]}"; do args+=(--vcf "$v"); done
  python3 "$SCRIPT_DIR/build_consensus_blacklist.py" \
    "${args[@]}" --fai "${REFERENCE}.fai" --min-support 2 \
    --union-output "$UNION_PLAIN" --consensus-output "$CONS_PLAIN"
  bgzip -f -@ "$THREADS" "$UNION_PLAIN"
  bgzip -f -@ "$THREADS" "$CONS_PLAIN"
  tabix -f -p vcf "${UNION_PLAIN}.gz"
  tabix -f -p vcf "${CONS_PLAIN}.gz"
  echo "Recommended conservative public WGS blacklist: ${CONS_PLAIN}.gz"
  echo "Broad flag-only union blacklist: ${UNION_PLAIN}.gz"
fi
