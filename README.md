# REWIRE RNA Editing Pipeline

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
│   ├── extract_high_confidence_candidates.py
│   └── make_wiki_figures.py
└── results/
    ├── CU5.15_EGFP_CC.control_corrected.PUF_candidates.tsv.gz
    └── CU5.15_EGFP_CC.high_confidence_candidates.tsv.gz
```

## Current result

| Category | Count |
|---|---:|
| Control-corrected candidate sites | 33,983 |
| High-confidence candidate sites | 2,955 |
| High-confidence sites with treated_edit_rate 0.10–0.20 | 2,446 |
| High-confidence sites with treated_edit_rate 0.20–0.50 | 503 |
| High-confidence sites with treated_edit_rate ≥0.50 | 6 |

## Pipeline summary

This workflow extracts transcriptome-wide C-to-U candidate editing sites from RNA-seq data for the REWIRE PUF-APOBEC system.

The workflow starts from SRA accessions, generates STAR-aligned BAM files, scans transcript-oriented C-equivalent sites in GENCODE exons, calculates editing rates, merges replicates at the count level, and performs treated-control correction.

Candidate editing sites are defined by elevated editing in the editor-treated group and low background editing in the matched mock/control group.
