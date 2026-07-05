# Cycle 1 — Build the RNA evidence branch

## Design

### Question

Can six RNA-seq libraries provide reproducible, transcript-oriented evidence for treatment-associated C-to-U editing?

### Experimental structure

| Condition | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

Each library was analysed independently. This preserved replicate support instead of pooling treated or control reads before calling.

### Evidence rules

The initial screening logic required:

```text
called in all three treated replicates
AND not called in any control replicate
AND consistent with transcript-level C-to-U editing
```

Transcript-level C-to-U is represented in genomic coordinates as:

```text
positive-strand transcript: genomic C>T
negative-strand transcript: genomic G>A
```

## Build

### Reference and paths

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
STAR_INDEX=/data/ydx/igem/STAR_index
REDITOOLS=/data/ydx/igem/REDItools2
```

The same GRCh38 FASTA was used for alignment, GATK preprocessing, VEP interpretation and subsequent catalogue validation.

### Stage 1 — FASTQ acquisition

```bash
python3 pipeline/scripts/rna/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest pipeline/config/samples.tsv \
  --threads 16
```

The downloader verifies that paired FASTQ files are present and non-empty before downstream processing.

### Stage 2 — STAR alignment

```bash
python3 pipeline/scripts/rna/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest pipeline/config/samples.tsv \
  --threads 50
```

STAR two-pass alignment was used so that splice junctions discovered in the first pass could inform the second pass.

### Stage 3 — GATK RNA preprocessing

```bash
python3 pipeline/scripts/rna/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest pipeline/config/samples.tsv \
  --java-options=-Xmx16g
```

The core operations were:

```text
MarkDuplicates
→ SplitNCigarReads
→ coordinate-sorted, indexed BAM files
```

Read-group information was checked because GATK fails when `@RG` records are absent.

### Stage 4 — REDItools2 calling

```bash
conda activate reditools2_py2

nohup bash pipeline/scripts/rna/run_reditools_all_samples.sh \
  "$PROJECT" "$REF" "$REDITOOLS" "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

Core REDItools2 settings:

| Setting | Role |
|---|---|
| `-S` | report positions containing substitutions |
| `-me 20` | require at least 20 edited reads for a reported call |
| mapping quality ≥20 | exclude poorly mapped reads |
| base quality ≥30 | exclude low-confidence bases |

The edited-read threshold was intentionally stringent. It favours strongly supported calls but may miss lower-frequency editing.

### Stage 5 — Union candidate coordinates

```bash
python3 pipeline/scripts/rna/reditools_union_to_vcf.py \
  --manifest pipeline/config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"
```

The union preserves every allele reported in at least one sample so that support can be compared across all six libraries.

### Stage 6 — Transcript orientation

```bash
python3 pipeline/scripts/rna/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache /path/to/vep_cache
```

VEP strand information was used to retain only events compatible with transcript-level C-to-U editing. Sites assigned to both transcript orientations were treated as ambiguous rather than forced into one class.

### Stage 7 — Replicate and control comparison

The integration code combined sample-level call status, editing rate, coverage, alternate-read count and VEP strand.

Current strict implementation:

```text
pipeline/scripts/rna/filter_c_to_u_and_compare.py
pipeline/scripts/rna/filter_calls.py
pipeline/scripts/rna/filter_utils.py
```

The frozen legacy result was produced with an earlier helper implementation, but the same treated-replicate and control-call logic was preserved.

## Test

### Per-sample strand-consistent calls

The frozen log reported the following values for five visible samples, with the sixth treated count recorded in the complete log:

```text
CU517_GC_T2    7,200
CU517_GC_T3    7,651
CU517_GC_C1    1,769
CU517_GC_C2    1,680
CU517_GC_C3    1,871
```

Treated samples contained substantially more strand-consistent calls than controls, supporting continued replicate-aware comparison rather than pooling.

### Evidence funnel before genomic filtering

```text
strand-consistent site matrix             9,930
called in all three treated replicates    4,778
treatment-specific before catalogue       3,349
```

The difference between 4,778 and 3,349 reflects candidates that were also called in at least one control.

## Learn

### Lesson 1 — A control non-call is not zero editing

REDItools2 `-me 20` means that a site can be absent from the call table even when the BAM contains lower-level alternate reads. Therefore:

```text
not called in control ≠ proven absence of editing
```

The strict workflow added direct candidate-site depth measurement in all six BAM files:

```bash
bash pipeline/scripts/rna/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  pipeline/config/samples.tsv
```

This script runs `samtools depth` with mapping quality 30 and base quality 20 over the union BED.

### Lesson 2 — The frozen legacy table has an evidence boundary

The completed legacy site matrix lacked the newer `all_replicates_depth_pass` column, and the original `candidate_depth/` directory was unavailable during final catalogue integration. We did not manufacture a replacement value.

Consequently, the frozen 3,349 pre-catalogue sites are described as treatment-associated screening candidates rather than fully depth-qualified events.

### Lesson 3 — Replicate-aware evidence should remain visible

The final tables retain individual treated coverage, alternate-read counts and editing rates. This allows later ranking by both editing magnitude and replicate consistency rather than relying on a pooled average alone.

## Output of this cycle

```text
CU5.17_EGFP_GC.site_matrix.tsv.gz
CU5.17_EGFP_GC.treated_consensus.tsv.gz
CU5.17_EGFP_GC.treatment_specific.tsv.gz
```

The treatment-specific table containing 3,349 candidates became the input to the genomic-catalogue cycles.
