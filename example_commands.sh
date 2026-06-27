#!/usr/bin/env bash
set -euo pipefail

# Example command set for the REWIRE RNA-seq to C-to-U editing site pipeline.
# Edit SRR IDs before running.

REF="/data/ydx/igem/GRCh38.primary_assembly.genome.fa"
GTF="/data/ydx/igem/gencode.v50.primary_assembly.annotation.gtf"
STAR_INDEX="/data/ydx/igem/STAR_index"
BAM_DIR="work_rewire_auto/bam"
PUF_TARGET="AACGUCUAUA"
THREADS=16

# 1. SRA to STAR BAM for treated replicates.
bash scripts/sra_to_star_bam_only.sh SRR_TREATED_1,SRR_TREATED_2 "$THREADS"

# 2. SRA to STAR BAM for mock/control replicates.
bash scripts/sra_to_star_bam_only.sh SRR_CONTROL_1,SRR_CONTROL_2 "$THREADS"

# 3. Generate treated all-C table by merging treated replicates at count level.
python3 scripts/all_genomic_c_editing_to_3files.py \
  --srr SRR_TREATED_1,SRR_TREATED_2 \
  --sample-id CU5.15_EGFP_CC \
  --editor CU5.15 \
  --puf-target "$PUF_TARGET" \
  --ref "$REF" \
  --gtf "$GTF" \
  --bam-dir "$BAM_DIR" \
  --out-prefix CU5.15_EGFP_CC \
  --min-coverage 20 \
  --min-baseq 20 \
  --edited-rate-cutoff 0.05 \
  --strong-rate-cutoff 0.5

# 4. Generate control all-C table by merging control replicates at count level.
python3 scripts/all_genomic_c_editing_to_3files.py \
  --srr SRR_CONTROL_1,SRR_CONTROL_2 \
  --sample-id control \
  --editor mock \
  --puf-target "$PUF_TARGET" \
  --ref "$REF" \
  --gtf "$GTF" \
  --bam-dir "$BAM_DIR" \
  --out-prefix control \
  --min-coverage 20 \
  --min-baseq 20 \
  --edited-rate-cutoff 0.05 \
  --strong-rate-cutoff 0.5

# 5. Combine treated and control into control-corrected candidate sites.
# This version avoids pandas/scipy and is easier to run on servers with restricted conda/pip.
python3 scripts/combine_treated_control_puf_candidates_nopandas.py \
  --treated CU5.15_EGFP_CC.all_C_sites.genomic.cov20.tsv.gz \
  --control control.all_C_sites.genomic.cov20.tsv.gz \
  --out CU5.15_EGFP_CC.control_corrected.PUF_candidates.tsv.gz \
  --min-treated-coverage 20 \
  --min-control-coverage 20 \
  --min-treated-rate 0.05 \
  --max-control-rate 0.02 \
  --min-delta 0.05 \
  --min-corrected-rate 0.05 \
  --ref "$REF" \
  --puf-target "$PUF_TARGET" \
  --flank 50 \
  --max-puf-mismatch 10
