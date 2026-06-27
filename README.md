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
│   └── make_wiki_figures.py
└── results/
    ├── CU5.15_EGFP_CC.control_corrected.PUF_candidates.tsv.gz
    └── CU5.15_EGFP_CC.high_confidence_candidates.tsv.gz
```

## Pipeline summary

This workflow extracts transcriptome-wide C-to-U candidate editing sites from RNA-seq data for the REWIRE PUF-APOBEC system.

The workflow starts from SRA accessions, generates STAR-aligned BAM files, scans transcript-oriented C-equivalent sites in GENCODE exons, calculates editing rates, merges replicates at the count level, and performs treated-control correction.

Candidate editing sites are defined by elevated editing in the editor-treated group and low background editing in the matched mock/control group.

## Main outputs

| Output | Description |
|---|---|
| `*.all_C_sites.genomic.cov20.tsv.gz` | All covered transcript-oriented C-equivalent sites |
| `*.edited_sites.genomic.cov20.rate005.tsv.gz` | Sites with editing rate ≥ 0.05 before control correction |
| `*.strong_sites.genomic.cov20.rate05.tsv.gz` | Sites with editing rate ≥ 0.5 before control correction |
| `*.control_corrected.PUF_candidates.tsv.gz` | Treated-control corrected candidate C-to-U editing sites |
| `*.high_confidence_candidates.tsv.gz` | Strictly filtered high-confidence candidate sites |

## Interpretation

The output sites are treated-control corrected C-to-U candidate editing sites. They should be used as candidate labels for downstream analysis, together with clean background C sites selected from the full `all_C_sites` table.
