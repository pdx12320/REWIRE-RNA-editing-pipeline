# RNA-seq-based C-to-U Editing Site Extraction Pipeline

## Purpose

The REWIRE PUF-APOBEC system edits RNA C to U. To identify transcriptome-wide candidate editing sites from RNA-seq data, we built a strand-aware site extraction pipeline using matched editor-treated and mock/control samples.

The pipeline extracts covered C-equivalent sites from exon regions, calculates site-level editing rates, merges biological replicates at the count level, and performs treated-control correction to reduce background signals.

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
GENCODE exon C-equivalent site scanning
↓
Base counting at each covered site
↓
Replicate count merging
↓
Treated-control correction
↓
Candidate C-to-U editing sites
↓
LAMAR training labels
```

## Strand-aware C-to-U signal definition

RNA-seq reads are aligned to the genomic reference. Therefore, transcript-level C-to-U editing appears as different genomic substitutions depending on transcript strand.

```text
+ strand transcript: genomic C→T
- strand transcript: genomic G→A
```

The negative-strand G→A signal is the genomic-coordinate representation of transcript-level C-to-U editing.

## Candidate editing site definition

A site is defined as a control-corrected candidate editing site if it satisfies:

```text
treated_total_count >= 20
control_total_count >= 20
treated_edit_rate >= 0.05
control_edit_rate <= 0.02
delta_edit_rate >= 0.05
corrected_edit_rate >= 0.05
```

High-confidence candidate sites can be selected using stricter thresholds:

```text
treated_edit_rate >= 0.10
corrected_edit_rate >= 0.10
control_edit_rate <= 0.02
```

If statistical testing is available, an additional threshold can be applied:

```text
Fisher q <= 0.05
```

## Directory structure

```text
rewire_editing_pipeline_wiki/
├── README.md
├── pipeline_method.md
├── results_summary.md
├── figures/
│   ├── pipeline_overview.png
│   ├── strand_definition.png
│   ├── treated_vs_control_scatter.png
│   ├── edit_rate_distribution.png
│   └── candidate_level_barplot.png
├── tables/
│   ├── parameter_table.tsv
│   ├── output_file_description.tsv
│   └── candidate_statistics.tsv
├── scripts/
│   ├── sra_to_star_bam_only.sh
│   ├── all_genomic_c_editing_to_3files.py
│   ├── combine_treated_control_puf_candidates.py
│   └── make_wiki_figures.py
└── results/
    ├── CU5.15_EGFP_CC.control_corrected.PUF_candidates.tsv.gz
    └── CU5.15_EGFP_CC.high_confidence_candidates.tsv.gz
```

## Output files

| File | Description | Use |
|---|---|---|
| `treated.all_C_sites.genomic.cov20.tsv.gz` | All covered C-equivalent sites in editor-treated samples | Candidate and background site pool |
| `control.all_C_sites.genomic.cov20.tsv.gz` | All covered C-equivalent sites in mock/control samples | Background correction |
| `*.control_corrected.PUF_candidates.tsv.gz` | Treated-control corrected candidate editing sites | Candidate positive set |
| `*.high_confidence_candidates.tsv.gz` | Strictly filtered candidate editing sites | High-confidence positives and result visualization |

## Model label construction

Positive samples can be selected from high-confidence candidate editing sites. Negative samples can be sampled from clean background C sites satisfying:

```text
treated_edit_rate < 0.02
control_edit_rate < 0.02
total_count >= 20
```

This converts transcriptome-wide RNA-seq editing signals into a candidate editing classification or ranking dataset.
