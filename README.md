# REWIRE Model 1

## Transcriptome-wide evidence for C-to-U RNA editing

REWIRE is designed to recruit a cytidine deaminase to a selected RNA target. Model 1 extends the analysis beyond the reporter and asks whether treatment-associated C-to-U signals can be distinguished from control background, insufficient coverage and plausible genomic variation.

![RNA-editing evidence generation pipeline](wiki/assets/figure1_model1_evidence_pipeline.svg)

The workflow combines two evidence branches:

- **RNA-seq:** three treated libraries and three controls are aligned, processed and screened for reproducible, transcript-oriented C-to-U signals.
- **Public WGS:** candidate HEK293T genome datasets are aligned and used to construct an external genomic-SNV blacklist.

The branches are integrated by exact `CHROM:POS:REF:ALT` matching. The final output is an auditable site-level table rather than an unqualified list of off-targets.

## Evidence required for a candidate

A conservative candidate must satisfy all of the following:

```text
called in all three treated replicates
not called in the three controls
covered by at least 20 reads in all six RNA-seq libraries
consistent with transcript-level C-to-U editing
absent from the selected exact-allele WGS blacklist
```

Public WGS is used as an external HEK293T genomic-variant catalogue, not as WGS matched to the exact experimental cell batch.

## Explore the project

| Resource | Contents |
|---|---|
| [iGEM-ready Model 1 page](wiki/README.md) | Scientific rationale, assumptions, workflow, validation, limitations and references |
| [Editable workflow figure](wiki/assets/figure1_model1_evidence_pipeline.svg) | Vector figure for the team wiki or presentation |
| [Pipeline implementation](pipeline/README.md) | Complete execution order and commands |
| [Expected outputs](pipeline/OUTPUTS.md) | RNA-seq, WGS and final evidence files |
| [Troubleshooting](pipeline/TROUBLESHOOTING.md) | Read-group, MPI, coverage and contig-name issues |

## Reproducible implementation

```text
pipeline/
├── config/       fixed RNA-seq and WGS manifests
├── env/          Conda environment definitions
└── scripts/
    ├── rna/      SRA, STAR, GATK, REDItools2, VEP and evidence integration
    └── wgs/      metadata checks, BWA alignment and bcftools SNV calling
```

Raw FASTQ, BAM, coverage, VCF and result files are excluded from version control. The repository contains the analysis logic required to regenerate them.
