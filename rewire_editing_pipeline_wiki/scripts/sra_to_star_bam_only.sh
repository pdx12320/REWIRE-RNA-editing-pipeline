#!/usr/bin/env bash
set -euo pipefail

SRR_LIST="${1:-}"
THREADS="${2:-16}"

if [[ -z "$SRR_LIST" ]]; then
  echo "Usage: bash sra_to_star_bam_only.sh SRR1,SRR2 16"
  exit 1
fi

STAR_INDEX="/data/ydx/igem/STAR_index"
WORKDIR="work_rewire_auto"
SRA_DIR="$WORKDIR/sra"
FASTQ_DIR="$WORKDIR/fastq"
BAM_DIR="$WORKDIR/bam"
LOG_DIR="$WORKDIR/logs"

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$BAM_DIR" "$LOG_DIR"
IFS=',' read -ra SRRS <<< "$SRR_LIST"

for SRR in "${SRRS[@]}"; do
  SRR="$(echo "$SRR" | xargs)"
  echo "========== $SRR =========="

  STAR_BAM="$BAM_DIR/${SRR}.Aligned.sortedByCoord.out.bam"
  STAR_BAI="${STAR_BAM}.bai"
  COMPAT_BAM="$BAM_DIR/${SRR}.split.bam"
  COMPAT_BAI="${COMPAT_BAM}.bai"

  if [[ ! -f "$STAR_BAM" ]]; then
    echo "[1/3] prefetch $SRR"
    if [[ ! -d "$SRA_DIR/$SRR" && ! -f "$SRA_DIR/$SRR.sra" ]]; then
      prefetch "$SRR" -O "$SRA_DIR"
    fi

    echo "[2/3] fasterq-dump $SRR"
    if [[ ! -f "$FASTQ_DIR/${SRR}_1.fastq.gz" && ! -f "$FASTQ_DIR/${SRR}.fastq.gz" ]]; then
      if [[ -f "$SRA_DIR/$SRR/$SRR.sra" ]]; then
        SRA_INPUT="$SRA_DIR/$SRR/$SRR.sra"
      elif [[ -f "$SRA_DIR/$SRR.sra" ]]; then
        SRA_INPUT="$SRA_DIR/$SRR.sra"
      else
        SRA_INPUT="$SRR"
      fi
      fasterq-dump "$SRA_INPUT" --split-files --threads "$THREADS" -O "$FASTQ_DIR"
      gzip -f "$FASTQ_DIR/${SRR}"*.fastq
    fi

    echo "[3/3] STAR align $SRR"
    if [[ -f "$FASTQ_DIR/${SRR}_1.fastq.gz" && -f "$FASTQ_DIR/${SRR}_2.fastq.gz" ]]; then
      STAR --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$FASTQ_DIR/${SRR}_1.fastq.gz" "$FASTQ_DIR/${SRR}_2.fastq.gz" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattrRGline ID:"$SRR" SM:"$SRR" PL:ILLUMINA LB:"$SRR" PU:"$SRR" \
        --outFileNamePrefix "$BAM_DIR/${SRR}." \
        > "$LOG_DIR/${SRR}.star.stdout.log" 2> "$LOG_DIR/${SRR}.star.stderr.log"
    elif [[ -f "$FASTQ_DIR/${SRR}.fastq.gz" ]]; then
      STAR --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$FASTQ_DIR/${SRR}.fastq.gz" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattrRGline ID:"$SRR" SM:"$SRR" PL:ILLUMINA LB:"$SRR" PU:"$SRR" \
        --outFileNamePrefix "$BAM_DIR/${SRR}." \
        > "$LOG_DIR/${SRR}.star.stdout.log" 2> "$LOG_DIR/${SRR}.star.stderr.log"
    else
      echo "Cannot find FASTQ for $SRR"
      exit 1
    fi
  fi

  if [[ ! -f "$STAR_BAI" ]]; then
    samtools index -@ "$THREADS" "$STAR_BAM"
  fi

  rm -f "$COMPAT_BAM" "$COMPAT_BAI"
  ln -s "$(basename "$STAR_BAM")" "$COMPAT_BAM"
  ln -s "$(basename "$STAR_BAI")" "$COMPAT_BAI"

done
