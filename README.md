# REWIRE Model 1: RNA-editing evidence pipeline

Model 1 converts treated/control RNA-seq and public HEK293T whole-genome sequencing (WGS) data into an auditable evidence matrix for candidate transcript-level C-to-U editing.

The repository is organised into two parts:

```text
wiki/       Copy-ready text and editable figures for the iGEM team wiki
pipeline/   Reproducible manifests, environments, scripts and technical notes
```

## Start here

| Goal | File |
|---|---|
| Copy Model 1 to the iGEM wiki | [`wiki/README.md`](wiki/README.md) |
| Download the editable workflow figure | [`wiki/assets/figure1_model1_evidence_pipeline.svg`](wiki/assets/figure1_model1_evidence_pipeline.svg) |
| Run or audit the analysis | [`pipeline/README.md`](pipeline/README.md) |
| Inspect expected outputs | [`pipeline/OUTPUTS.md`](pipeline/OUTPUTS.md) |
| Resolve known deployment issues | [`pipeline/TROUBLESHOOTING.md`](pipeline/TROUBLESHOOTING.md) |

## Evidence design

```text
RNA-seq branch
3 treated + 3 controls
→ STAR
→ GATK RNA preprocessing
→ coverage-aware REDItools2
→ VEP transcript orientation
→ all-sample depth and control subtraction

WGS branch
3 public HEK293T WGS runs
→ metadata check
→ BWA-MEM2/BWA-MEM
→ GATK MarkDuplicates
→ bcftools SNV calling
→ merged or 2-of-3 exact-allele blacklist

Integration
replicate evidence + control evidence + strand evidence + WGS overlap
→ treatment-associated C-to-U candidates
```

The public WGS data are an external HEK293T genomic-variant catalogue. They are not WGS matched to the exact CU5.17 experimental cell batch.

## Repository tree

```text
.
├── README.md
├── LICENSE
├── wiki/
│   ├── README.md
│   └── assets/
│       └── figure1_model1_evidence_pipeline.svg
└── pipeline/
    ├── README.md
    ├── OUTPUTS.md
    ├── TROUBLESHOOTING.md
    ├── config/
    │   ├── samples.tsv
    │   └── wgs_runs.tsv
    ├── env/
    │   ├── reditools2_py2.yml
    │   └── wgs_pipeline.yml
    └── scripts/
        ├── rna/
        └── wgs/
```

Large FASTQ, BAM, coverage, VCF and result files are intentionally excluded from GitHub.
