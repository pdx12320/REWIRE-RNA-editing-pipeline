# REWIRE RNA-seq to C-to-U Editing Site Pipeline

This repository contains a Wiki-ready dry-lab package for extracting transcriptome-wide C-to-U candidate editing sites from RNA-seq/SRA data and converting them into AI-ready labels for downstream LAMAR-style modeling.

## Purpose

The pipeline was designed for the PUF-APOBEC / REWIRE system. Starting from raw SRA files, it produces strand-aware C-to-U editing tables, performs treated-control correction, and exports candidate editing sites that can be used as positive labels for model training.

The key idea is not to directly treat all apparent edited sites in the treated sample as true PUF editing. Instead, each covered C-equivalent site is compared between the editor-treated group and a matched mock/control group.

## Repository structure

```text
REWIRE-RNA-editing-pipeline/
├── README.md
├── pipeline_method.md
├── results_summary_template.md
├── example_commands.sh
├── figures/
│   ├── pipeline_overview.svg
│   └── strand_definition.svg
└── scripts/
    ├── sra_to_star_bam_only.sh
    ├── all_genomic_c_editing_to_3files.py
    └── combine_treated_control_puf_candidates_nopandas.py
```

## Workflow

```text
SRA files
↓
FASTQ extraction
↓
STAR alignment to GRCh38
↓
Sorted BAM files
↓
Scan GENCODE exon C-equivalent sites
↓
Count strand-aware C-to-U signals
↓
Merge treated replicates and mock/control replicates at count level
↓
Treated-control correction
↓
Control-corrected candidate editing sites
↓
LAMAR training labels
```

## Editing signal definition

Because RNA-seq reads are aligned to the genomic reference, C-to-U editing appears differently depending on transcript strand.

```text
+ strand transcript: genomic C→T
- strand transcript: genomic G→A
```

The negative-strand G→A signal does not mean that the editor edits G. It is the genomic-coordinate representation of transcript-level C-to-U editing.

## Main output files

| File | Meaning | Downstream use |
|---|---|---|
| `treated.all_C_sites.genomic.cov20.tsv.gz` | all covered transcript-oriented C-equivalent sites in treated samples | background pool and model input |
| `control.all_C_sites.genomic.cov20.tsv.gz` | all covered transcript-oriented C-equivalent sites in mock/control samples | background correction |
| `*.control_corrected.PUF_candidates.tsv.gz` | treated-control corrected candidate editing sites | positive candidate set |
| `high_confidence_candidates.tsv.gz` | stricter candidate subset | Wiki figures and model positives |

## Recommended candidate definition

A site is considered a control-corrected candidate if it satisfies:

```text
treated_total_count >= 20
control_total_count >= 20
treated_edit_rate >= 0.05
control_edit_rate <= 0.02
delta_edit_rate >= 0.05
corrected_edit_rate >= 0.05
```

For high-confidence positives, use stricter filters, for example:

```text
treated_edit_rate >= 0.10
corrected_edit_rate >= 0.10
control_edit_rate <= 0.02
```

If Fisher exact test results are available, add:

```text
fisher_q <= 0.05
```

## How this connects to LAMAR

Positive samples should be selected from control-corrected or high-confidence candidate editing sites. Negative samples should be sampled from clean background C sites, for example:

```text
treated_edit_rate < 0.02
control_edit_rate < 0.02
total_count >= 20
```

This makes the modeling task a candidate editing classifier/ranking model, rather than a direct precise editing-efficiency regression task.

## Important interpretation note

The output sites should be described as **control-corrected C-to-U candidate editing sites**, not as fully confirmed PUF-edited sites. PUF target similarity is best used as an auxiliary annotation or downstream enrichment analysis, not as a hard filter.
