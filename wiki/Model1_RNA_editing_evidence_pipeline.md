# Model 1 — RNA-editing evidence pipeline

## Overview

Our REWIRE system is designed to recruit a cytidine deaminase to a selected RNA sequence. Editing at the intended reporter site demonstrates that the construct is active, but it does not establish transcriptome-wide specificity. To evaluate broader editing evidence, we developed a reproducible RNA-seq workflow that compares three editor-treated libraries with three matched controls.

The purpose of Model 1 is not to label every RNA-seq mismatch as an off-target. Instead, it builds an evidence chain. A candidate must first be supported by quality-filtered reads, then interpreted in transcript orientation, reproduced across treated replicates, evaluated at the same coordinate in controls, and optionally checked against matched genomic variation.

## Biological question

Model 1 addresses three questions:

1. Which substitutions are reproducibly detected after editor expression?
2. Which substitutions are consistent with transcript-level C-to-U editing?
3. Which candidates remain after control, depth, strand, and genomic-variant filtering?

## Experimental design

We analyzed six paired-end RNA-seq libraries:

| Condition | Replicate | Sample ID | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

Replicates are central to the design. A mismatch detected in only one library may arise from sequencing error, mapping ambiguity, stochastic low-level noise, or sample-specific variation. Requiring independent support across treated replicates makes the final evidence more conservative and interpretable.

## Workflow

```text
SRA paired-end RNA-seq
→ FASTQ conversion
→ STAR two-pass alignment to GRCh38
→ read-group validation
→ GATK MarkDuplicates
→ GATK SplitNCigarReads
→ coverage-aware parallel REDItools2 calling
→ union substitution VCF
→ VEP transcript-strand annotation
→ transcript-oriented C-to-U interpretation
→ candidate-site depth confirmation in all six samples
→ treated/control replicate comparison
→ optional matched HEK293T WGS filtering
```

## Step 1 — Download and convert sequencing data

The six SRA accessions are downloaded with the NCBI SRA Toolkit. `fasterq-dump` converts each accession into paired FASTQ files, which are compressed before alignment.

```bash
python3 scripts/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest config/samples.tsv \
  --threads 16
```

This step creates one read-1 and one read-2 FASTQ file for every sample.

## Step 2 — Align reads to GRCh38

Paired-end reads are aligned to the GRCh38 primary assembly with STAR in two-pass mode. Two-pass alignment allows splice junctions detected in the first pass to inform the final alignment. Coordinate-sorted BAM files are produced for downstream analysis.

```bash
python3 scripts/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest config/samples.tsv \
  --threads 50
```

Read-group fields are written during alignment so that downstream GATK tools can identify the sample, library, platform, and read group.

## Step 3 — Preprocess RNA alignments with GATK

RNA-seq alignments require additional preparation before site-level mismatch analysis.

```bash
python3 scripts/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest config/samples.tsv \
  --java-options=-Xmx16g
```

The preprocessing stage performs three operations:

1. **Read-group validation or repair.** Missing read groups are added before duplicate processing.
2. **MarkDuplicates.** PCR duplicate reads are flagged and duplicate metrics are recorded.
3. **SplitNCigarReads.** Spliced RNA alignments are transformed into a representation compatible with downstream site-level analysis.

## Step 4 — Build the REDItools2 coverage map

The MPI implementation of REDItools2 requires a position-level coverage map before site calling. This map is used to divide the genome into intervals with more balanced computational workloads.

```bash
bash scripts/generate_reditools_coverage_limited.sh \
  SAMPLE.splitncigarreads.bam \
  SAMPLE_COVERAGE_DIRECTORY \
  GRCh38.primary_assembly.genome.fa.fai \
  8
```

The final argument limits the number of simultaneous `samtools depth` processes. This avoids launching one disk-intensive job for every reference contig and reduces I/O contention on a shared server.

GRCh38 supplementary contigs such as `GL000194.1` and `KI270750.1` are valid reference sequences. Their `.1` and `.2` suffixes are part of the contig identifiers and must be preserved.

## Step 5 — Call candidate substitutions with REDItools2

Each sample is processed separately, while genomic intervals within that sample are distributed across MPI workers.

```bash
conda activate reditools2_py2

nohup bash scripts/run_reditools_all_samples.sh \
  "$PROJECT" \
  "$REF" \
  "$REDITOOLS" \
  "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

The final three numbers specify:

```text
30  MPI processes for REDItools2
8   simultaneous coverage jobs
8   compression threads
```

REDItools2 is run with:

```text
-S       report edited sites only
-me 20   require at least 20 edited reads
```

The `-me 20` setting is an edited-read threshold. It is not equivalent to requiring total depth of 20. Because it is stringent, it prioritizes robust signals but can miss low-frequency editing events.

For each reported position, REDItools2 records:

```text
chromosome
position
reference base
strand field
quality-filtered coverage
mean base quality
A/C/G/T read counts
observed substitution type
estimated substitution frequency
```

## Step 6 — Combine substitutions into a union VCF

The six REDItools2 tables are combined into one union set of observed substitutions. A BED file containing the union candidate coordinates is also generated.

```bash
mkdir -p "$PROJECT/vcf"

python3 scripts/reditools_union_to_vcf.py \
  --manifest config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"
```

The union representation ensures that every candidate coordinate can be evaluated consistently across all six samples.

## Step 7 — Annotate transcript strand with VEP

A genomic substitution cannot be interpreted as transcript-level C-to-U editing without knowing transcript orientation. The union VCF is therefore annotated with Ensembl VEP in offline GRCh38 mode.

```bash
python3 scripts/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache "$VEP_CACHE"
```

The key rule is:

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

A genomic G→A call on a negative-strand transcript is the reverse-complement representation of transcript-level C-to-U editing. It is not interpreted as biochemical editing of G.

Coordinates with contradictory strand support across overlapping transcripts are treated as ambiguous and excluded rather than assigned an arbitrary orientation.

## Step 8 — Measure candidate-site depth in every replicate

Absence from a REDItools2 call table is not automatically evidence of absence. A site may simply be uncovered or insufficiently sequenced in that sample. We therefore measure depth independently at every union candidate coordinate in all six BAM files.

```bash
bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  config/samples.tsv
```

The depth query uses:

```text
base quality ≥30
mapping quality ≥20
```

This step distinguishes an adequately sequenced negative call from a missing observation.

## Step 9 — Construct the evidence matrix

The final comparison integrates REDItools2 calls, VEP strand, all-replicate depth, treated/control status, and optional WGS overlap.

```bash
python3 scripts/filter_c_to_u_and_compare.py \
  --manifest config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

When matched HEK293T WGS data are available, the following option is added:

```text
--wgs-vcf /path/to/HEK293T.filtered.vcf.gz
```

## Conservative default definition

A treatment-specific candidate must satisfy all of the following:

```text
reported in all three treated replicates
reported in no control replicate
candidate-site depth ≥20 in all six samples
consistent with transcript-level C-to-U orientation
absent from the matched WGS variant set, when available
```

The full per-sample matrix is retained rather than exporting only a final yes/no label. Each row stores call status, candidate-site depth, REDItools2 depth, edited-read count, and editing frequency for every replicate. This makes the filtering logic transparent and allows alternative thresholds to be tested without repeating the upstream pipeline.

## Quality-control logic

The workflow includes several safeguards:

- BAM files are checked with `samtools quickcheck`.
- Read-group metadata are validated before GATK duplicate processing.
- Coverage maps and REDItools2 temporary intervals are stored separately for each sample.
- Parallel outputs are sorted according to the reference FASTA index before merging.
- Contig version suffixes such as `.1` and `.2` are preserved.
- Final compressed tables are indexed with tabix.
- Candidate depth is measured independently from REDItools2 call status.
- Treated and control evidence is retained at the replicate level.

## What Model 1 contributes

Model 1 is the evidence-generation layer of the REWIRE dry-lab framework. Its role is not merely to run an editing caller. It converts raw sequencing data into a traceable, replicate-aware candidate matrix suitable for biological interpretation and downstream modeling.

This evidence layer is important because a downstream machine-learning model can only be as reliable as its labels. By separating alignment, site discovery, strand interpretation, control comparison, and coverage confirmation, Model 1 provides a more defensible source of candidate positives and background sites for later prediction tasks.

## Interpretation and limitations

The final output should be described as a set of **computational RNA-editing candidates**, not automatically confirmed biological off-targets. RNA-seq mismatches may also arise from:

- genomic variants;
- alignment ambiguity;
- sequencing error;
- endogenous RNA modification;
- library-specific artifacts;
- low-complexity or repetitive sequence.

Matched controls, replicate consistency, depth confirmation, and WGS filtering reduce these alternatives but cannot eliminate them completely. High-priority candidates require independent validation, such as targeted amplicon sequencing, independent RNA-seq, Sanger sequencing for strong signals, or another orthogonal assay.

## Results status

The complete workflow and code are available in the repository. Numerical result files, final candidate counts, and result figures are intentionally omitted until all six samples complete the same analysis and pass integrity checks.
