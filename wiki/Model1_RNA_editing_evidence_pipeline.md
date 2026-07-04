# Model 1 — RNA-editing Evidence Pipeline

## From raw sequencing reads to an auditable C-to-U evidence matrix

REWIRE recruits a cytidine deaminase to a selected RNA sequence. Reporter editing demonstrates activity at the intended target, but specificity requires a transcriptome-wide analysis. Model 1 therefore asks:

> **Which C-to-U signals are reproducible in treated samples, absent or minimal in controls, supported by sufficient read depth, consistent with transcript orientation, and not readily explained by genomic variation?**

The complete copy-ready Methods text is available in [`Model1_Methods_copy_ready.md`](Model1_Methods_copy_ready.md). Source code, environments and technical notes are available in the [GitHub repository](https://github.com/pdx12320/REWIRE-RNA-editing-pipeline).

---

## Figure 1 — RNA-editing evidence generation pipeline

![Figure 1. RNA-editing evidence generation pipeline](assets/model1_pipeline_with_wgs.svg)

**Figure 1. RNA-editing evidence generation pipeline.** The RNA-seq branch identifies quality-supported substitutions and evaluates transcript orientation, all-sample depth and treated/control reproducibility. The WGS branch aligns three public HEK293T whole-genome sequencing runs, calls genomic SNVs and produces either a merged call set or an exact-allele two-of-three consensus blacklist. Both branches converge in a site-level evidence matrix. Public WGS is used as an external blacklist rather than as WGS matched to the exact experimental cell batch.

---

## Experimental design

Three editor-treated RNA-seq libraries were compared with three control libraries.

| Condition | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

Replicates test whether a signal is reproducible rather than library-specific. Controls provide background context, but only when the same coordinate has enough sequencing depth to be informative.

Three public HEK293T WGS runs are configured as an external genomic-variant resource:

```text
SRR37832939
SRR37832940
SRR37832941
```

Before processing, their ENA metadata are checked to determine whether they represent multiple runs from one BioSample or independent public HEK293T genomes.

---

# RNA-seq branch

## Stage 1 — Download and organize RNA-seq data

The treated/control assignment is fixed in `config/samples.tsv`. SRA Toolkit downloads each accession and converts it into paired FASTQ files.

```bash
python3 scripts/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest config/samples.tsv \
  --threads 16
```

**Checkpoint:** every sample must produce non-empty read-1 and read-2 FASTQ files.

## Stage 2 — STAR two-pass alignment

Paired-end reads are aligned to the GRCh38 primary assembly with STAR. Two-pass alignment first discovers splice junctions and then uses them during final alignment. Coordinate-sorted BAM files and read-group fields are written for each sample.

```bash
python3 scripts/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest config/samples.tsv \
  --threads 50
```

Core STAR configuration:

```bash
STAR \
  --genomeDir "$STAR_INDEX" \
  --runThreadN 50 \
  --readFilesIn sample_1.fastq.gz sample_2.fastq.gz \
  --readFilesCommand gunzip -c \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattrRGline ID:sample SM:sample LB:sample PL:ILLUMINA PU:sample
```

Read groups identify the sample, library and platform. They do not normalize sequencing depth; they allow GATK and Picard to process each sample and library correctly.

## Stage 3 — GATK RNA preprocessing

GATK MarkDuplicates flags PCR/optical duplicates and records duplicate metrics. SplitNCigarReads then processes reads spanning splice junctions into exon-aligned segments suitable for mismatch analysis.

```bash
python3 scripts/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest config/samples.tsv \
  --java-options=-Xmx16g
```

**Checkpoint:** every final BAM is validated with `samtools quickcheck` and has an index.

## Stage 4 — Coverage map for REDItools2 parallelization

Parallel REDItools2 uses a per-position coverage map to divide the reference into intervals with more balanced computational loads. This coverage map is a scheduling input, not the final editing result.

```bash
bash scripts/generate_reditools_coverage_limited.sh \
  sample.splitncigarreads.bam \
  sample_coverage_directory \
  GRCh38.primary_assembly.genome.fa.fai \
  8
```

Coverage generation is limited to eight concurrent `samtools depth` processes to reduce disk contention. REDItools2 itself uses 30 MPI processes.

## Stage 5 — REDItools2 substitution discovery

Each sample is analyzed independently, while genomic intervals within each sample are distributed across MPI workers.

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

Core REDItools2 command:

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

| Parameter | Function |
|---|---|
| `-S` | output positions with observed substitutions |
| `-me 20` | require at least 20 edited reads at a reported site |
| default `-q 20` | minimum mapping quality |
| default `-bq 30` | minimum base quality |
| `-G`, `-D` | complete and per-contig coverage inputs |
| `-np 30` | 30 MPI processes |

The `-me 20` value is an edited-read threshold, not a total-depth threshold. It prioritizes strongly supported events but may miss low-frequency activity.

### Engineering fix — preserve complete GRCh38 contig identifiers

GRCh38 includes supplementary contigs such as `GL000194.1`, `GL000205.2` and `KI270750.1`. Their `.1` and `.2` suffixes are part of the reference identifier. The original temporary-file parser removed everything after the first dot and failed during final sorting. We replaced it with:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

This removes only the final `.gz` extension and preserves the exact contig name.

## Stage 6 — VEP transcript-strand interpretation

REDItools2 reports substitutions in genomic coordinates, while C-to-U editing occurs in transcripts. Substitutions from all six samples are combined into one union VCF and BED file, then annotated with VEP.

```bash
python3 scripts/reditools_union_to_vcf.py \
  --manifest config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"

python3 scripts/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache "$VEP_CACHE"
```

Interpretation rule:

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

The negative-strand G→A representation is the reverse complement of transcript-level C-to-U editing; it does not imply biochemical G editing. Coordinates assigned to transcripts on both orientations are treated as ambiguous.

## Stage 7 — Independent depth in all six samples

A missing REDItools2 call is not automatically evidence of no editing. The candidate coordinate may simply have insufficient coverage. We therefore query every union candidate in every final BAM using base quality at least 30 and mapping quality at least 20.

```bash
bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  config/samples.tsv
```

This separates:

```text
not called despite sufficient depth
from
not observed because sequencing depth was insufficient
```

---

# Public HEK293T WGS branch

## Stage 8 — Check WGS metadata

The WGS pipeline first queries ENA metadata. Consecutive SRR identifiers alone do not prove that the runs belong to one sample.

```bash
bash scripts/wgs/00_check_sra_metadata.sh \
  config/wgs_runs.tsv \
  wgs_sra_metadata.tsv
```

The script checks:

```text
library_strategy = WGS
library_layout = PAIRED
scientific_name = Homo sapiens
sample_accession = same or different across runs
```

If all runs share one BioSample, their BAM files are merged and called together. If they represent different BioSamples, each run is called separately and exact alleles supported by at least two of three call sets form the conservative blacklist.

## Stage 9 — WGS alignment and SNV calling

The complete WGS workflow is:

```bash
conda env create -f environment/wgs_pipeline.yml
conda activate rewire_wgs

REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
WGS_OUT=/data/ydx/igem/HEK293T_public_WGS_3runs

nohup bash scripts/wgs/run_3run_wgs_pipeline.sh \
  --runs config/wgs_runs.tsv \
  --reference "$REF" \
  --outdir "$WGS_OUT" \
  --mode auto \
  --threads 32 \
  --min-dp 10 \
  --min-alt 3 \
  --min-vaf 0.05 \
  --min-qual 20 \
  > "$WGS_OUT.pipeline.log" 2>&1 &
```

The branch performs:

```text
SRA Toolkit download and conversion
→ BWA-MEM2/BWA-MEM alignment to the same GRCh38 FASTA
→ coordinate sorting and indexing
→ GATK MarkDuplicates
→ bcftools mpileup and variant calling
→ normalization and SNP-only filtering
→ merged call set or exact-allele 2-of-3 consensus blacklist
```

Default single-run WGS thresholds:

```text
depth ≥10
alternate reads ≥3
alternate-allele fraction ≥0.05
QUAL ≥20
FILTER is PASS or unset
```

Recommended outputs:

```text
same BioSample:
  HEK293T_3runs.filtered.SNV.vcf.gz

different BioSamples:
  HEK293T_3runs.consensus2of3.SNV.vcf.gz   conservative exclusion list
  HEK293T_3runs.union.SNV.vcf.gz           broad annotation list
```

These public runs are an external HEK293T genomic-variant catalogue, not matched WGS from the exact CU5.17 experimental cell batch.

---

# Evidence integration

## Stage 10 — Control subtraction and final evidence matrix

The final table integrates:

- transcript-oriented C-to-U status;
- treated and control REDItools2 calls;
- candidate-site depth in all six BAMs;
- edited-read counts and editing fractions;
- overlap with the selected exact-allele WGS blacklist.

```bash
python3 scripts/filter_c_to_u_and_compare.py \
  --manifest config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --wgs-vcf "$WGS_OUT/vcf/HEK293T_3runs.consensus2of3.SNV.vcf.gz" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

The conservative default definition requires:

```text
called in all three treated replicates
AND called in no control replicate
AND candidate-site depth ≥20 in all six samples
AND transcript orientation consistent with C-to-U editing
AND exact CHROM:POS:REF:ALT absent from the selected WGS blacklist
```

WGS comparison is exact-allele based. A different alternate base at the same coordinate does not automatically remove the RNA candidate.

---

# What Model 1 produces

Model 1 retains an auditable evidence matrix rather than only a final yes/no label. For every candidate and sample, the matrix records:

```text
coordinate
reference and alternate allele
transcript orientation
REDItools2 call status
independent candidate-site depth
edited-read count
editing fraction
treated/control replicate support
WGS blacklist status
```

This structure supports alternative thresholds without rerunning alignment or REDItools2.

---

# Connection to the wet lab

Model 1 does not replace experimental validation. It prioritizes sites for targeted amplicon sequencing, independent RNA-seq or another orthogonal assay. It can also provide evidence-derived labels for downstream sequence models: reproducible candidates become positive examples, while sufficiently covered background positions can form a controlled negative set.

---

# Limitations

RNA-seq mismatches can arise from genomic variants, alignment ambiguity, sequencing artifacts, endogenous RNA modification, repetitive regions and batch effects. Treated/control comparison, independent depth and public HEK293T WGS filtering reduce these alternatives but do not eliminate them.

Because the WGS data are public and not from the exact experimental cell batch, retained sites should be described as:

> **treatment-associated RNA-editing candidates filtered against an external HEK293T genomic-variant catalogue**

not as definitively SNV-free editing events. High-priority candidates still require orthogonal validation.

---

# Our contribution

Model 1 provides:

- a fixed three-treated/three-control RNA-seq design;
- STAR and GATK RNA preprocessing;
- coverage-aware MPI REDItools2 discovery;
- preservation of versioned GRCh38 contig identifiers;
- VEP-based transcript-oriented C-to-U interpretation;
- independent depth assessment in all six samples;
- treated/control replicate filtering;
- a reproducible three-run public HEK293T WGS workflow;
- merged or two-of-three consensus genomic-SNV blacklists;
- an auditable site-level evidence matrix;
- reusable code, environments and copy-ready Wiki text.

Numerical result files are intentionally excluded until all samples complete the same workflow and pass the same integrity checks.
