# REWIRE RNA Editing Pipeline (Model 1)

This repository documents the **RNA-editing evidence pipeline** used in the REWIRE project. The goal is to derive high-confidence C-to-U (and G-to-A on reverse strand) RNA editing candidates from RNA-seq data using a reproducible multi-step workflow.

---

## Overview

The pipeline converts raw SRA samples into curated RNA editing candidate tables through the following stages:

1. SRA → FASTQ → BAM (STAR alignment)
2. SplitNCigarReads (GATK RNA preprocessing)
3. Coverage construction (REDItools2 requirement)
4. REDItools2 parallel scanning (MPI)
5. Interval merging and table construction
6. Treated vs control filtering (future stage)

---

## Input Data

- RNA-seq samples (paired treated/control)
- Reference genome: GRCh38 primary assembly
- STAR index (prebuilt)
- Sample manifest (`config/samples.tsv`)

Example SRA sets:
- Treated: SRR27885768, SRR27885766, SRR27885765
- Control: SRR27885767, SRR27885764, SRR27885763

---

## Core Pipeline Script

Main execution entry:

```bash
/data/ydx/igem/run_cu517_egfp_gc_paper_pipeline_v2.sh
```

Or updated modular version:

```bash
scripts/run_reditools_all_samples.sh
```

---

## Step 1: Alignment (STAR)

```bash
STAR \
  --runThreadN 16 \
  --genomeDir STAR_index \
  --readFilesIn sample.fastq \
  --outSAMtype BAM SortedByCoordinate
```

Output:
- sorted BAM
- indexed BAM

---

## Step 2: GATK SplitNCigarReads

Used to prepare RNA-seq alignments for variant calling.

```bash
gatk SplitNCigarReads \
  -R reference.fa \
  -I input.bam \
  -O split.bam
```

---

## Step 3: Coverage Generation

Coverage is required by REDItools2 to accelerate base-level scanning.

```bash
samtools depth split.bam > coverage.txt
```

Notes:
- includes all chromosomes (chr1–chrY, GL*, KI*)
- must preserve contig version numbers (.1, .2)

---

## Step 4: REDItools2 Parallel Scan

MPI-based scanning of RNA mismatches.

```bash
mpirun -np 30 python2 parallel_reditools.py \
  -f split.bam \
  -r reference.fa \
  -G coverage_file \
  -D coverage_dir \
  -t temp_dir \
  -Z reference.fa.fai \
  -S -me 20
```

Key parameters:
- `-S`: strand-aware mode
- `-me 20`: minimum coverage threshold
- `-np`: MPI workers

---

## Known Implementation Fixes

### 1. Contig parsing bug (critical)

Fixed issue:
- GL/KI contig names were truncated incorrectly

Fix applied:
```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

### 2. Missing read group / MPI errors

Resolved by ensuring BAM compatibility and correct MPI slot allocation.

---

## Step 5: Merge REDItools2 Output

All interval chunks are merged into a single table:

- Region
- Position
- Base counts
- Editing frequency

Output format:
```
Region Position Reference Strand Coverage ... Frequency
```

---

## Step 6: Treated vs Control Filtering

Final candidate definition:

- Treated editing rate ≥ threshold
- Control background rate ≤ threshold
- Minimum coverage ≥ 20
- Optional Fisher test filtering

---

## Output Structure

```
reditools/
├── coverage/
├── tmp/
├── tables/
│   ├── *_T1.txt.gz
│   ├── *_T2.txt.gz
│   └── *_high_confidence.tsv.gz (future)
└── logs/
```

---

## Notes

- GL/KI contigs are retained with version suffixes (.1, .2)
- Pipeline is MPI-parallel and cluster-scale
- Coverage stage is required for REDItools2 interval splitting
- All intermediate files are preserved for reproducibility

---

## Status

- T1 (REDItools2 scan): completed
- T2–T6: in progress
- Final candidate set: pending
