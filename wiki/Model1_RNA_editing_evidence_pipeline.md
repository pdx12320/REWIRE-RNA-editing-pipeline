# Model 1 — RNA-editing Evidence Pipeline

## From RNA-seq reads to evidence for REWIRE activity

### Why did we build Model 1?

Our REWIRE system recruits a cytidine deaminase to a selected RNA sequence. Editing at the designed reporter can show that the construct is active, but a reporter alone cannot answer an equally important question:

> **Does the editor also produce reproducible C-to-U signals elsewhere in the transcriptome?**

We therefore built Model 1 as an evidence-generation pipeline. Instead of treating every RNA-seq mismatch as an off-target, Model 1 asks whether a candidate is supported by high-quality reads, reproduced across treated samples, adequately observed in controls, consistent with transcript orientation, and separable from possible genomic variation.

Model 1 is also the foundation for later sequence models. A machine-learning model is only useful when its training labels are credible; this pipeline creates the evidence from which those labels can be built.

---

## Our experimental design

We used three editor-treated RNA-seq libraries and three matched controls.

| Condition | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

The replicate structure gives Model 1 two kinds of evidence:

- **reproducibility**, because a candidate should not depend on one library;
- **background context**, because a missing or weak control signal changes how a treated call is interpreted.

<!-- FIGURE 1 PLACEHOLDER
Six sample cards: T1, T2, T3 on the left and C1, C2, C3 on the right.
Use one visual style for treated samples and another for controls.
-->

---

## The Model 1 workflow

```text
Six RNA-seq libraries
        ↓
SRA download and FASTQ conversion
        ↓
STAR alignment to GRCh38
        ↓
GATK RNA preprocessing
        ↓
Coverage-aware parallel REDItools2 analysis
        ↓
Per-sample substitution tables
        ↓
VEP transcript-strand annotation
        ↓
Transcript-oriented C-to-U interpretation
        ↓
Depth confirmation in all six samples
        ↓
Treated/control evidence matrix
        ↓
Candidates for validation and downstream modeling
```

<!-- FIGURE 2 PLACEHOLDER
Full-width horizontal workflow.
Recommended labels: Input, Map, Prepare, Call, Orient, Compare, Validate.
-->

---

# Stage 1 — Organize the samples

The relationship between sample name, condition, replicate number, and SRA accession is stored in one manifest:

```text
config/samples.tsv
```

```tsv
sample	group	replicate	srr
CU517_GC_T1	treated	1	SRR27885768
CU517_GC_T2	treated	2	SRR27885766
CU517_GC_T3	treated	3	SRR27885765
CU517_GC_C1	control	1	SRR27885767
CU517_GC_C2	control	2	SRR27885764
CU517_GC_C3	control	3	SRR27885763
```

Keeping this information in a machine-readable file reduces manual sample handling and lets every downstream script use the same treated/control assignment.

The reads are downloaded and converted into paired compressed FASTQ files:

```bash
python3 scripts/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest config/samples.tsv \
  --threads 16
```

**Checkpoint:** every sample must produce one non-empty read-1 file and one non-empty read-2 file.

---

# Stage 2 — Map the reads

We aligned paired-end reads to the GRCh38 primary assembly with STAR in two-pass mode.

```bash
python3 scripts/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest config/samples.tsv \
  --threads 50
```

The STAR command writes coordinate-sorted BAM files and includes read-group information for each sample.

```bash
STAR \
  --genomeDir "$STAR_INDEX" \
  --runThreadN 50 \
  --readFilesIn sample_1.fastq.gz sample_2.fastq.gz \
  --readFilesCommand gunzip -c \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattrRGline ID:sample SM:sample PL:ILLUMINA LB:sample PU:sample
```

### Why read groups became part of our design

During development, missing read-group metadata caused GATK to fail. We therefore turned read-group presence into a quality-control checkpoint rather than treating it as an invisible software detail.

```bash
samtools view -H sample.Aligned.sortedByCoord.out.bam | grep '^@RG'
```

**Checkpoint:** the BAM must be readable, indexed, and linked to the correct sample through its `@RG` line.

---

# Stage 3 — Prepare RNA alignments

RNA-seq reads often span exon junctions, so we processed each alignment before site calling.

```bash
python3 scripts/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest config/samples.tsv \
  --java-options=-Xmx16g
```

The preprocessing stage performs three operations:

1. repair read groups when necessary;
2. mark PCR duplicates and record duplicate metrics;
3. run `SplitNCigarReads` on spliced RNA alignments.

Representative command:

```bash
gatk SplitNCigarReads \
  -R GRCh38.primary_assembly.genome.fa \
  -I markduplicates.bam \
  -O splitncigarreads.bam
```

We mark duplicates rather than silently deleting them. This keeps the decision visible in the BAM flags and preserves the original read evidence for later inspection.

**Checkpoint:** every SplitNCigarReads BAM is tested with `samtools quickcheck` before entering REDItools2.

---

# Stage 4 — Build the REDItools2 coverage map

The parallel REDItools2 implementation needs a per-position coverage map. This file is not the final editing result. It is used to divide the genome into computational intervals with more balanced workloads.

```bash
bash scripts/generate_reditools_coverage_limited.sh \
  sample.splitncigarreads.bam \
  sample_coverage_directory \
  GRCh38.primary_assembly.genome.fa.fai \
  8
```

The helper script calculates coverage separately for each reference contig:

```bash
samtools depth -r chromosome sample.splitncigarreads.bam \
  > sample_coverage_directory/chromosome
```

### Why we limited the number of coverage jobs

The original helper can launch one background process for every reference contig. On shared storage, this can create many disk-intensive jobs at the same time. We separated the two forms of parallelism:

```text
coverage generation: 8 concurrent depth jobs
REDItools2 analysis: 30 MPI processes
```

This retains parallel analysis while reducing avoidable I/O contention.

---

# Stage 5 — Call candidate substitutions

Each sample is analyzed separately with the MPI implementation of REDItools2.

```bash
nohup bash scripts/run_reditools_all_samples.sh \
  "$PROJECT" \
  "$REF" \
  "$REDITOOLS" \
  "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

The main REDItools2 call is:

```bash
mpirun -np 30 python2 parallel_reditools.py \
  -f sample.splitncigarreads.bam \
  -r GRCh38.primary_assembly.genome.fa \
  -G sample.splitncigarreads.cov \
  -D sample_coverage_directory/ \
  -t sample_temp_directory/ \
  -Z GRCh38.primary_assembly.genome.fa.fai \
  -S \
  -me 20
```

| Parameter | Role in Model 1 |
|---|---|
| `-S` | strict mode: only positions with an observed edit are written |
| `-me 20` | require at least 20 editing events at a reported position |
| default `-q 20` | remove reads below mapping quality 20 |
| default `-bq 30` | remove bases below base quality 30 |
| `-np 30` | distribute intervals across 30 MPI processes |
| `-G`, `-D` | supply the complete and per-contig coverage data |

The edited-read threshold is deliberately stringent. It prioritizes strongly supported sites but may miss genuine low-frequency events.

For every reported position, REDItools2 records the coordinate, reference base, quality-filtered coverage, A/C/G/T counts, observed substitution, and estimated frequency.

---

# Engineering challenge — Preserving GRCh38 contig names

GRCh38 contains supplementary contigs with versioned names such as:

```text
GL000194.1
GL000205.2
KI270750.1
```

The `.1` and `.2` suffixes are part of the reference identifiers.

During our first full run, REDItools2 completed the interval computation but failed during output sorting:

```text
ValueError: 'chrGL000009' is not in list
```

The original filename parser removed everything after the first dot. This changed:

```text
GL000009.2#100#500.gz
```

into:

```text
GL000009
```

which no longer matched the reference FASTA index.

We replaced the parser with:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

This removes only the final `.gz` extension and preserves the complete contig name:

```python
["GL000009.2", "100", "500"]
```

<!-- FIGURE 3 PLACEHOLDER
Before/after debugging diagram:
GL000009.2 → GL000009 (wrong)
GL000009.2 → GL000009.2 (correct)
-->

### What we learned

A pipeline can complete hours of biological computation and still fail during file integration. Reference naming and merge logic therefore became explicit parts of our validation process.

---

# Stage 6 — Interpret C-to-U in transcript orientation

The BAM and REDItools2 tables use genomic coordinates, but RNA editing must be interpreted in transcript orientation.

```text
positive-strand transcript:
transcript C-to-U = genomic C-to-T
```

```text
negative-strand transcript:
transcript C-to-U = genomic G-to-A
```

A genomic G-to-A call on a negative-strand transcript does not mean that the editor biochemically edits G. It is the reverse-complement representation of transcript-level C-to-U editing.

We combine substitutions from all six sample tables into a union VCF and candidate BED:

```bash
python3 scripts/reditools_union_to_vcf.py \
  --manifest config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"
```

We then annotate transcript strand with offline VEP:

```bash
python3 scripts/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache "$VEP_CACHE"
```

Coordinates with contradictory strand annotations across overlapping transcripts are treated as ambiguous instead of being assigned arbitrarily.

<!-- FIGURE 4 PLACEHOLDER
Two-panel strand figure:
+ transcript: C→T
− transcript: G→A
Both panels point to “transcript-level C→U”.
-->

---

# Stage 7 — Ask whether every sample had enough evidence

A site absent from a control call table is not automatically a true negative. It may simply have too little coverage.

To separate these cases, we query every union candidate position in every BAM using:

```text
base quality ≥30
mapping quality ≥20
```

```bash
bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  config/samples.tsv
```

This produces an independent depth value even when REDItools2 did not report an edit.

### Why a missing call can mean different things

```text
Control A: no edit call, depth = 125
Interpretation: the position was observed but did not pass the edit-call rule
```

```text
Control B: no edit call, depth = 0
Interpretation: the sample provides no evidence about this position
```

Model 1 keeps these two situations separate.

<!-- FIGURE 5 PLACEHOLDER
Three evidence cards:
1. treated call + covered control negative;
2. treated call + uncovered control;
3. call in only one treated replicate.
-->

---

# Stage 8 — Build the six-sample evidence matrix

The final comparison combines:

- REDItools2 call status;
- VEP transcript orientation;
- candidate-site depth in every sample;
- treated/control assignment;
- optional overlap with a matched WGS variant set.

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

Our conservative default evidence rule is:

```text
called in all three treated replicates
AND called in no control replicate
AND candidate-site depth ≥20 in all six samples
AND consistent with transcript-level C-to-U editing
AND absent from the optional matched WGS variant set
```

The final matrix stores more than a pass/fail result. For each sample, it retains call status, independent depth, REDItools2 depth, edited-read count, and editing frequency. This lets us test alternative thresholds later without rerunning alignment and site calling.

<!-- FIGURE 6 PLACEHOLDER
Filtering funnel:
all substitutions → strand-consistent C-to-U → treated reproducibility
→ control exclusion → all-sample depth → optional WGS exclusion.
Do not add numerical counts until the complete analysis is available.
-->

---

# How Model 1 supports the wet lab

Model 1 is not intended to replace experimental validation. It reduces the search space and records why each site was prioritized.

## 1. Prioritizing validation targets

Candidates with reproducible treated support and informative control coverage can be ranked for targeted amplicon sequencing, independent RNA-seq, or another orthogonal assay.

## 2. Evaluating specificity

The intended reporter site and endogenous candidates can be summarized using the same evidence types: depth, edited-read support, frequency, replicate consistency, and control status.

## 3. Creating labels for downstream models

High-confidence positives and well-covered background positions can become training examples for later sequence models. Model 1 provides the evidence layer; the later model learns patterns from that evidence.

<!-- FIGURE 7 PLACEHOLDER
Wet Lab–Dry Lab feedback loop:
RNA-seq → Model 1 → ranked candidates → targeted validation → improved labels.
-->

---

# Design–Build–Test–Learn

## Design

We designed a six-sample workflow based on replicate consistency, matched controls, transcript orientation, and all-sample depth.

## Build

We combined SRA Toolkit, STAR, GATK, samtools, MPI REDItools2, VEP, and custom Python filters into a modular pipeline.

## Test

We tested BAM readability, read-group presence, coverage generation, interval completion, compressed output integrity, reference ordering, and tabix indexing.

## Learn

Three lessons changed the final design:

1. a missing call is not evidence without enough coverage;
2. transcript orientation changes the genomic representation of C-to-U editing;
3. exact reference-contig names must survive every temporary filename and merge operation.

---

# Limitations

Model 1 produces computational candidates, not automatically confirmed off-targets.

RNA-seq mismatches can also arise from:

- genomic variants;
- alignment ambiguity;
- sequencing artifacts;
- endogenous RNA modification;
- repetitive or low-complexity regions;
- batch-specific effects.

Matched HEK293T WGS filtering is important for stronger conclusions. Without matched WGS, retained positions should be described as **RNA-derived candidate editing sites**.

The strict edited-read threshold may also miss low-frequency activity. A sensitivity analysis can later examine lower thresholds together with stronger artifact controls.

High-priority sites still require orthogonal validation.

---

# Our contribution

Model 1 provides:

- a fixed six-sample treated/control manifest;
- a reproducible RNA-seq-to-evidence workflow;
- limited-concurrency coverage generation for shared servers;
- a fix for REDItools2 versioned-contig parsing;
- transcript-oriented C-to-U interpretation;
- independent depth confirmation in every sample;
- an auditable treated/control evidence matrix;
- reusable code for future REWIRE constructs.

Complete commands are available in [`Code_and_commands.md`](Code_and_commands.md). Final numerical results and figures will be added only after all six samples finish the same workflow and pass the same checks.
