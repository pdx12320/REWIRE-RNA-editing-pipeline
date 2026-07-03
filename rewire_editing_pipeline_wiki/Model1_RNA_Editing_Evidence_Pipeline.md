# Model 1: RNA-editing Evidence Pipeline

## From RNA-seq reads to evidence for REWIRE activity

### The question behind Model 1

Our REWIRE system is designed to recruit a cytidine deaminase to a selected RNA sequence. A reporter experiment can show that the construct works at the intended target, but it cannot answer a second question that is equally important:

> **Does the editor also introduce reproducible C-to-U changes at other RNA sites?**

To address this question, we developed Model 1: a transcriptome-wide RNA-editing evidence pipeline. The pipeline does not begin with a prediction model. It begins with experimental RNA-seq data and builds a traceable chain of evidence from raw sequencing reads to replicate-supported candidate sites.

This distinction matters. A downstream machine-learning model can only be as reliable as the labels used to train it. Model 1 therefore serves as the evidence-generation layer for the rest of our Dry Lab work.

---

## Our design

We compared three editor-treated RNA-seq libraries with three matched control libraries.

| Condition | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

The three-replicate design allows us to distinguish reproducible signals from mismatches that appear in only one library. Controls allow us to identify background substitutions that are unrelated to the introduced editor.

### Workflow at a glance

```text
Raw RNA-seq data
        ↓
SRA download and FASTQ conversion
        ↓
STAR alignment to GRCh38
        ↓
GATK RNA preprocessing
        ↓
Coverage-aware REDItools2 analysis
        ↓
Per-sample mismatch tables
        ↓
VEP transcript-strand annotation
        ↓
Transcript-oriented C-to-U interpretation
        ↓
Depth confirmation in all six samples
        ↓
Treated/control replicate comparison
        ↓
Candidate sites for validation and downstream modeling
```

<!-- WIKI FIGURE PLACEHOLDER
Figure 1. A horizontal workflow diagram using one icon per stage.
Recommended labels: RNA-seq, Alignment, RNA preprocessing, Site calling,
Strand annotation, Replicate comparison, Candidate evidence.
-->

---

# 1. Building a consistent input dataset

We stored the relationship between sample name, experimental group, replicate number, and SRA accession in a single manifest file:

```text
config/samples.tsv
```

This prevents sample names from being entered manually at different stages and reduces the risk of mixing treated and control libraries.

```tsv
sample	group	replicate	srr
CU517_GC_T1	treated	1	SRR27885768
CU517_GC_T2	treated	2	SRR27885766
CU517_GC_T3	treated	3	SRR27885765
CU517_GC_C1	control	1	SRR27885767
CU517_GC_C2	control	2	SRR27885764
CU517_GC_C3	control	3	SRR27885763
```

The raw sequencing data are downloaded with the SRA Toolkit and converted into paired FASTQ files.

```bash
python3 scripts/download_sra_fastq.py \
  --project /data/ydx/igem/CU5.17_EGFP_GC_paper \
  --manifest config/samples.tsv \
  --threads 16
```

### Checkpoint

Before alignment, each sample must have two non-empty compressed FASTQ files:

```text
sample_1.fastq.gz
sample_2.fastq.gz
```

---

# 2. Mapping reads to the human reference genome

We aligned paired-end reads to the GRCh38 primary assembly with STAR in two-pass mode. Two-pass alignment allows STAR to use splice junctions detected in the first pass when performing the second alignment pass.

```bash
python3 scripts/run_star_alignment.py \
  --project /data/ydx/igem/CU5.17_EGFP_GC_paper \
  --star-index /data/ydx/igem/STAR_index \
  --manifest config/samples.tsv \
  --threads 50
```

The underlying STAR command writes a coordinate-sorted BAM file and includes sample read-group tags:

```bash
STAR \
  --genomeDir "$STAR_INDEX" \
  --runThreadN 50 \
  --readFilesIn sample_1.fastq.gz sample_2.fastq.gz \
  --readFilesCommand gunzip -c \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattrRGline \
      ID:sample SM:sample PL:ILLUMINA LB:sample PU:sample
```

### Why read groups matter

GATK uses read-group information to identify the sample, library, platform, and sequencing unit associated with each read. During development, missing read groups caused a GATK failure. We therefore made read-group validation an explicit checkpoint rather than assuming the BAM header was complete.

### Checkpoint

```bash
samtools quickcheck -v sample.Aligned.sortedByCoord.out.bam
samtools view -H sample.Aligned.sortedByCoord.out.bam | grep '^@RG'
```

A valid BAM must pass `samtools quickcheck`, have a BAM index, and contain the expected `@RG` line.

---

# 3. Preparing RNA alignments for site-level analysis

RNA-seq reads frequently cross exon junctions. We therefore processed each STAR BAM with GATK before running REDItools2.

The three operations are:

1. repair read groups when they are missing;
2. mark PCR duplicates;
3. run `SplitNCigarReads` to process spliced RNA alignments.

```bash
python3 scripts/run_gatk_preprocessing.py \
  --project /data/ydx/igem/CU5.17_EGFP_GC_paper \
  --reference /data/ydx/igem/GRCh38.primary_assembly.genome.fa \
  --manifest config/samples.tsv \
  --java-options=-Xmx16g
```

Representative commands:

```bash
gatk AddOrReplaceReadGroups \
  -I star.bam \
  -O readgroups.bam \
  --RGID sample \
  --RGLB sample \
  --RGPL ILLUMINA \
  --RGPU sample \
  --RGSM sample
```

```bash
gatk MarkDuplicates \
  -I readgroups.bam \
  -O markduplicates.bam \
  -M markduplicates_metrics.txt \
  --CREATE_INDEX true
```

```bash
gatk SplitNCigarReads \
  -R GRCh38.primary_assembly.genome.fa \
  -I markduplicates.bam \
  -O splitncigarreads.bam
```

### What we keep

Duplicates are marked rather than silently discarded. This preserves the original alignments and makes the duplicate decision visible in the BAM flags and metrics file.

### Checkpoint

```bash
samtools quickcheck -v sample.splitncigarreads.bam
```

Only indexed, readable SplitNCigarReads BAM files are passed to REDItools2.

---

# 4. Creating a coverage map for parallel analysis

Parallel REDItools2 does not use the coverage map as the final editing result. Instead, it uses per-position read depth to divide the genome into intervals with more balanced computational workloads.

The reference FASTA contains the main chromosomes as well as supplementary contigs with names such as:

```text
GL000194.1
GL000205.2
KI270750.1
```

The version suffix is part of the contig identifier and must be preserved.

We calculate coverage separately for each contig and limit the number of simultaneous `samtools depth` processes to avoid overwhelming shared storage.

```bash
bash scripts/generate_reditools_coverage_limited.sh \
  sample.splitncigarreads.bam \
  sample_coverage_directory \
  GRCh38.primary_assembly.genome.fa.fai \
  8
```

The helper script performs the equivalent of:

```bash
samtools depth -r chromosome sample.splitncigarreads.bam \
  > sample_coverage_directory/chromosome
```

and then concatenates all contig files in FASTA-index order.

### Why we limited concurrency

The original coverage helper can start one background process for every contig. On a shared server, this can create many simultaneous disk-intensive jobs. We separated coverage parallelism from REDItools2 MPI parallelism:

```text
Coverage generation: 8 concurrent depth jobs
REDItools2 analysis: 30 MPI processes
```

This gives the alignment data time to be read efficiently without removing the computational benefit of MPI during the main analysis.

---

# 5. Calling candidate mismatches with REDItools2

We run the parallel REDItools2 implementation separately for each of the six samples.

```bash
bash scripts/run_reditools_all_samples.sh \
  /data/ydx/igem/CU5.17_EGFP_GC_paper \
  /data/ydx/igem/GRCh38.primary_assembly.genome.fa \
  /data/ydx/igem/REDItools2 \
  /path/to/reditools2_py2/bin/python \
  30 8 8
```

The central MPI command is:

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

### Parameters used in Model 1

| Parameter | Meaning in our workflow |
|---|---|
| `-S` | strict mode; only positions containing an observed edit are written |
| `-me 20` | require at least 20 editing events at a reported position |
| default `-q 20` | discard reads below mapping quality 20 |
| default `-bq 30` | exclude bases below base quality 30 |
| `-np 30` | use 30 MPI processes |
| `-G` and `-D` | provide the complete coverage file and per-contig coverage directory |

The `-me 20` threshold is not the same as requiring total depth 20. It requires edited-read support at the REDItools2 discovery stage and is therefore deliberately stringent.

### Per-sample output

The merged table contains:

```text
Region
Position
Reference
Strand
Coverage-q30
MeanQ
BaseCount[A,C,G,T]
AllSubs
Frequency
five genomic-DNA columns
```

The genomic-DNA columns remain empty because RNA and DNA are not analyzed together inside REDItools2 in this workflow. Genomic-variant filtering is handled later with an optional matched WGS VCF.

---

# 6. A bug we found and fixed

During the first complete run, REDItools2 finished the interval calculations but failed while sorting temporary output files:

```text
ValueError: 'chrGL000009' is not in list
```

The cause was not a missing chromosome. The original code removed everything after the first dot in a filename. As a result:

```text
GL000009.2#start#end.gz
```

was incorrectly converted to:

```text
GL000009
```

but the GRCh38 FASTA index contains the full contig name `GL000009.2`.

We changed the parser from removing everything after the first dot to removing only the final `.gz` extension:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

Now the parser correctly produces:

```python
["GL000009.2", "start", "end"]
```

This fix preserves every version suffix, including `.1`, `.2`, and any future suffix included in the reference.

### What this taught us

A pipeline can complete hours of computation and still fail during file integration. Reproducibility therefore requires checking not only biological thresholds but also reference naming, file naming, and merge logic.

<!-- WIKI CARD PLACEHOLDER
Title: Debugging contribution
Text: We identified and repaired a REDItools2 filename-parsing bug that affected
versioned GRCh38 supplementary contigs. The fix preserves the exact FASTA contig ID.
-->

---

# 7. Converting genomic mismatches into transcript-level C-to-U evidence

RNA editing is interpreted in transcript orientation, but the BAM and REDItools2 coordinates are genomic.

For a transcript on the positive strand:

```text
transcript C-to-U → genomic C-to-T
```

For a transcript on the negative strand:

```text
transcript C-to-U → genomic G-to-A
```

The negative-strand G-to-A signal does not indicate biochemical editing of G. It is the reverse-complement representation of C-to-U editing in the RNA transcript.

We combine the substitutions reported across all six samples into a union VCF:

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
  --cache /path/to/vep_cache
```

We retain only substitutions consistent with transcript-oriented C-to-U editing. Coordinates receiving contradictory positive- and negative-strand transcript annotations are treated as ambiguous and are not assigned arbitrarily.

<!-- WIKI FIGURE PLACEHOLDER
Figure 2. Two-panel strand diagram.
Left: + transcript, genomic C→T.
Right: − transcript, genomic G→A.
Caption: Both signals represent transcript-level C-to-U conversion.
-->

---

# 8. Checking evidence in every replicate

A site absent from a REDItools2 output table is not automatically a true negative. It may simply have insufficient coverage in that sample.

To avoid this mistake, we query every union candidate position in every SplitNCigarReads BAM using:

```text
minimum base quality = 30
minimum mapping quality = 20
```

```bash
bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  config/samples.tsv
```

This creates an independent depth value even when REDItools2 did not report an editing event.

### Why this is important

Consider two control samples:

```text
Control A: no edit call, depth = 125
Control B: no edit call, depth = 0
```

Control A provides evidence that the position was sequenced without meeting the editing-call threshold. Control B provides no evidence about the position because it was not covered.

Model 1 keeps these cases separate.

---

# 9. Building the treated/control evidence matrix

The final filtering step brings all six samples together.

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

Our conservative default candidate rule is:

```text
called in all 3 treated replicates
AND called in 0 control replicates
AND depth ≥20 in every treated and control replicate
AND consistent with transcript-level C-to-U editing
AND absent from the optional matched WGS variant set
```

The matrix preserves more than a pass/fail label. For each sample it stores:

- whether REDItools2 called the site;
- independent candidate-site depth;
- REDItools2 quality-filtered depth;
- edited-read count;
- estimated editing frequency.

This means the team can later test less stringent or more stringent criteria without rerunning the computationally expensive alignment and REDItools2 stages.

---

# 10. How Model 1 supports the wet lab

Model 1 was not designed as an isolated bioinformatics exercise. It supports experimental decisions in three ways.

## Prioritizing validation targets

Sites supported by all treated replicates, adequately covered in controls, and absent from controls can be ranked for targeted amplicon sequencing or another independent assay.

## Distinguishing on-target activity from background

The pipeline keeps the reporter/on-target site and endogenous candidate sites in the same evidence framework: genomic coordinate, edited-read support, depth, replicate consistency, and control comparison.

## Creating labels for downstream modeling

High-confidence positives and well-covered background sites can be used to construct a downstream sequence-learning dataset. In this project, Model 1 provides the evidence layer, while later models learn sequence patterns from the labels generated here.

<!-- WIKI FIGURE PLACEHOLDER
Figure 3. Wet Lab–Dry Lab feedback loop.
Wet lab RNA-seq → Model 1 evidence → ranked sites → targeted validation → improved labels.
-->

---

# 11. Design–Build–Test–Learn

## Design

We designed a six-sample treated/control workflow that required replicate consistency, transcript-strand interpretation, and independent depth confirmation.

## Build

We implemented a modular pipeline around SRA Toolkit, STAR, GATK, samtools, MPI REDItools2, VEP, and custom Python filters.

## Test

At each stage we checked file readability, BAM indexing, read-group presence, coverage-file generation, interval completion, compressed-table integrity, and tabix indexing.

## Learn

We learned that three details strongly affect the reliability of RNA-editing evidence:

1. absence of a call is not evidence without adequate depth;
2. transcript strand changes the genomic representation of C-to-U editing;
3. exact reference-contig names must be preserved through parallel file generation and merging.

These lessons were fed back into the final pipeline design.

---

# 12. Limitations

Model 1 reduces technical and biological alternatives, but it does not turn every retained site into a confirmed off-target.

RNA-seq mismatches may still arise from:

- genomic variants;
- alignment ambiguity;
- sequencing artifacts;
- endogenous RNA modifications;
- unmodeled batch effects;
- low-complexity or repetitive regions.

Matched HEK293T WGS filtering is therefore important for stronger claims. Without matched WGS, retained sites should be described as **RNA-derived candidate editing sites**, not definitive editing-only events.

The strict `-me 20` discovery setting also favors strongly supported sites and may miss genuine low-frequency editing. A later sensitivity analysis can lower this threshold while applying stronger artifact controls.

Finally, high-priority candidates require orthogonal validation. Suitable options include targeted amplicon sequencing, independent RNA-seq, or Sanger sequencing for high-frequency events.

---

# 13. Our contribution

Model 1 contributes more than a one-time analysis:

- a fixed six-sample manifest;
- a reproducible RNA-seq-to-evidence workflow;
- a coverage-limited REDItools2 parallelization strategy;
- a fix for versioned GL/KI contig parsing;
- explicit distinction between genomic mismatches and transcript-level C-to-U editing;
- independent all-replicate depth confirmation;
- an auditable treated/control evidence matrix;
- reusable code for future REWIRE constructs and datasets.

The complete implementation is available in the repository `scripts/` directory. Numerical results and final figures will be added only after all six samples complete the same workflow and pass the same integrity checks.
