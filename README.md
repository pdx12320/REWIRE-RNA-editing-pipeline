# REWIRE Model 1

## Transcriptome-wide evidence for C-to-U RNA editing

REWIRE recruits a programmable PUF–APOBEC editor to a selected RNA target. Model 1 asks whether treatment-associated C-to-U signals can be separated from control background, insufficient coverage and plausible genomic variation. It combines replicate-aware RNA-seq evidence with an assembly-harmonized 293T variant catalogue from the [HEK293 Genome Project](https://hek293genome.org/v2/data.php).

![RNA-editing evidence generation pipeline](wiki/assets/figure1_model1_evidence_pipeline.svg)

The workflow has two evidence branches:

- **RNA evidence:** three treated libraries and three controls are aligned, processed and screened independently for reproducible, transcript-oriented C-to-U signals.
- **293T genomic catalogue:** the database-released `293T_CG` call set is converted from NCBI build 36/hg18 to GRCh38, validated against the same reference FASTA and used as an exact-allele exclusion flag.

The branches are integrated by exact `CHROM:POS:REF:ALT` matching. Catalogue-overlapping sites remain in the complete site matrix but are excluded from the high-confidence treatment-specific set.

## Evidence required for a high-confidence candidate

```text
called in all three treated replicates
not called in the three controls
covered by at least 20 reads in all six RNA-seq libraries
consistent with transcript-level C-to-U editing
absent from the selected 293T genomic catalogue
```

The catalogue is external evidence, not whole-genome sequencing matched to the exact experimental cell batch. The final sites should therefore be described as **treatment-associated RNA-editing candidates**, not definitively SNV-free off-targets.

## Why the catalogue workflow replaced public-SRA WGS reconstruction

Three candidate public WGS BioSamples were evaluated, but only 19.2–26.3% of reads mapped to GRCh38 and the two-of-three consensus contained only 118 variants. This coverage was insufficient for a genome-wide exclusion resource. The final workflow therefore uses the database-released `293T_CG` catalogue and makes the assembly conversion explicit and reproducible.

## Explore the project

| Resource | Contents |
|---|---|
| [iGEM-ready Model 1 page](wiki/README.md) | Question, design decision, method, validation, results, limitations and references |
| [Editable workflow figure](wiki/assets/figure1_model1_evidence_pipeline.svg) | Vector figure for the team wiki or presentation |
| [Pipeline implementation](pipeline/README.md) | Complete execution order and commands |
| [Catalogue provenance](pipeline/CATALOGUE_PROVENANCE.md) | Source, assembly conversion, QC counts and interpretation boundary |
| [Expected outputs](pipeline/OUTPUTS.md) | RNA evidence, catalogue and integrated result files |
| [Troubleshooting](pipeline/TROUBLESHOOTING.md) | REDItools2, liftover, reference and indexing issues |

## Reproducible implementation

```text
pipeline/
├── config/                   fixed RNA-seq manifest
├── env/                      Conda environment definitions
└── scripts/
    ├── rna/                  SRA, STAR, GATK, REDItools2, VEP and evidence integration
    └── catalogue/            hg18-to-GRCh38 conversion of the 293T_CG catalogue
```

Raw FASTQ, BAM, coverage, VCF and result files are excluded from version control. The repository contains the analysis logic and provenance needed to regenerate the final evidence tables.
