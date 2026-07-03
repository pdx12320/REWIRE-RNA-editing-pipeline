# Model 1 — RNA-editing evidence pipeline

This repository contains the reproducible computational workflow used to generate RNA-editing evidence for the REWIRE project. It converts six paired-end RNA-seq libraries into an auditable matrix of transcript-oriented C-to-U candidate sites.

The repository contains the full workflow, source code, installation notes, troubleshooting records, and text that can be adapted directly to an iGEM Wiki. Numerical result tables are intentionally not committed while the six-sample analysis is still running.

## Scientific question

Reporter editing confirms that an engineered editor is active, but it does not establish transcriptome-wide specificity. Model 1 therefore asks:

1. Which RNA substitutions are reproducibly detected in editor-treated samples?
2. Which of those substitutions are consistent with transcript-level C-to-U editing?
3. Which candidates remain after matched-control, sequencing-depth, strand, and optional genomic-variant filtering?

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
→ depth confirmation in all six replicates
→ treated/control comparison
→ optional matched HEK293T WGS filtering
```

## Dataset

| Group | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

The machine-readable manifest is stored in `config/samples.tsv`.

## Repository layout

```text
config/
  samples.tsv                         fixed sample manifest

docs/
  INSTALLATION.md                     software and reference setup
  RUNBOOK.md                          operational notes
  STEP_BY_STEP.md                     numbered workflow
  PIPELINE.md                         methodological design
  OUTPUTS.md                          output schemas
  TROUBLESHOOTING.md                  deployment errors and fixes
  LIMITATIONS.md                      interpretation boundaries

scripts/
  download_sra_fastq.py               SRA download and FASTQ conversion
  run_star_alignment.py               STAR two-pass alignment
  run_gatk_preprocessing.py           GATK RNA preprocessing
  generate_reditools_coverage_limited.sh
  run_reditools_all_samples.sh        six-sample REDItools2 runner
  reditools_union_to_vcf.py           union VCF and candidate BED
  run_vep_annotation.py               transcript-strand annotation
  build_candidate_depth_tables.sh     all-replicate depth confirmation
  filter_utils.py
  filter_calls.py
  filter_c_to_u_and_compare.py        final evidence matrix

wiki/
  Model1_RNA_editing_evidence_pipeline.md
  Code_and_commands.md
  Methods.md
  Filtering.md
  Interpretation.md
  RESULTS_TEMPLATE.md

results/
  placeholders only; numerical results are not committed yet
```

## Main methodological decisions

### Strand-aware definition

A biochemical C-to-U event has two genomic representations:

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

The negative-strand G→A call is the reverse-complement representation of transcript-level C-to-U editing. It is not interpreted as biochemical G editing.

### Independent depth confirmation

A missing REDItools2 call is not automatically evidence of absence. The candidate coordinate is queried independently in every BAM using base quality 30 and mapping quality 20. This distinguishes an adequately sequenced negative sample from an uncovered site.

### Conservative default evidence rule

The default treatment-specific rule requires:

```text
called in all three treated replicates
called in no control replicate
candidate-site depth ≥20 in every replicate
strand consistent with transcript-level C-to-U editing
absent from an optional matched HEK293T WGS variant set
```

These thresholds are configurable and the full per-sample matrix is retained so alternative definitions can be evaluated without repeating alignment and site calling.

## REDItools2 implementation notes

Parallel REDItools2 requires both a combined coverage file and per-contig coverage files. The repository uses a limited-concurrency coverage generator to reduce disk contention on shared storage.

GRCh38 contigs such as `GL000194.1` and `KI270750.1` include version suffixes as part of their identifiers. The REDItools2 filename parser was patched so only the final `.gz` extension is removed; `.1`, `.2`, and other contig-version suffixes are preserved.

## Results status

The methods and executable workflow are complete. Numerical result files, final site counts, and result figures will be added only after all six samples finish the same workflow and pass integrity checks.

## iGEM Wiki material

The main copy-ready page is:

- `wiki/Model1_RNA_editing_evidence_pipeline.md`

Exact code blocks are provided in:

- `wiki/Code_and_commands.md`

The original repository state is preserved on branch `backup-before-paper-pipeline-20260704`.
