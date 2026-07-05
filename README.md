# REWIRE Model 1

## Transcriptome-wide evidence for C-to-U RNA editing

REWIRE recruits a programmable PUF–APOBEC editor to a selected RNA target. Model 1 asks whether treatment-associated C-to-U signals can be separated from control-call background and plausible genomic variation. It combines replicate-aware RNA-seq evidence with an assembly-harmonized 293T variant catalogue from the [HEK293 Genome Project](https://hek293genome.org/v2/data.php).

![RNA-editing evidence generation pipeline](wiki/assets/figure1_model1_evidence_pipeline.svg)

The workflow has two evidence branches:

- **RNA evidence:** three treated libraries and three controls are aligned, processed and screened independently for reproducible, transcript-oriented C-to-U signals.
- **293T genomic catalogue:** the database-released `293T_CG` call set is converted from NCBI build 36/hg18 to GRCh38, validated against the same reference FASTA and used as an exact-allele exclusion flag.

The branches are integrated by exact `CHROM:POS:REF:ALT` matching. Catalogue-overlapping sites remain auditable but are excluded from the retained screening set.

## Frozen Model 1 results

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact 293T catalogue overlaps | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

The final 3,333 sites are suitable for Lamar inference and candidate ranking after strand-oriented sequence extraction. They are not yet a complete supervised training dataset because missing control calls must not be treated as zero editing, and the frozen table does not contain independent depth/base-count evidence for those control non-calls.

## Implemented screening definition

```text
called in all three treated replicates
not called in the three controls under the original REDItools2 filter
consistent with transcript-level C-to-U editing
absent from the selected 293T genomic catalogue
```

A stricter high-confidence definition would additionally require direct candidate-site depth and base counts in all six BAM files. The catalogue is external evidence, not whole-genome sequencing matched to the exact experimental cell batch. The final sites should therefore be described as **catalogue-filtered treatment-associated screening candidates**, not definitively SNV-free off-targets.

## Why the catalogue workflow replaced public-SRA WGS reconstruction

Three candidate public WGS BioSamples were evaluated, but only 19.2–26.3% of reads mapped to GRCh38 and the two-of-three consensus contained only 118 variants. This coverage was insufficient for a genome-wide exclusion resource. The final workflow therefore uses the database-released `293T_CG` catalogue and makes the assembly conversion explicit and reproducible.

## Explore the project

| Resource | Contents |
|---|---|
| [iGEM-ready Model 1 page](wiki/README.md) | Question, design decision, method, frozen results, evidence boundary, limitations and references |
| [Frozen result summary](results/README.md) | The 9,930 → 4,778 → 3,349 → 3,333 evidence funnel |
| [Lamar handoff](model2/README.md) | Sequence orientation, provisional labels and inference/training boundary |
| [Editable workflow figure](wiki/assets/figure1_model1_evidence_pipeline.svg) | Vector figure for the team wiki or presentation |
| [Pipeline implementation](pipeline/README.md) | Complete execution order, legacy-output compatibility and commands |
| [Catalogue provenance](pipeline/CATALOGUE_PROVENANCE.md) | Source, assembly conversion, QC counts and interpretation boundary |
| [Expected outputs](pipeline/OUTPUTS.md) | RNA evidence, catalogue, integration and Lamar handoff files |
| [Troubleshooting](pipeline/TROUBLESHOOTING.md) | REDItools2, legacy-output, liftover, reference and indexing issues |

## Reproducible implementation

```text
pipeline/
├── config/                   fixed RNA-seq manifest
├── env/                      Conda environment definitions
└── scripts/
    ├── rna/                  SRA, STAR, GATK, REDItools2 and VEP
    ├── catalogue/            hg18-to-GRCh38 conversion and legacy-table filtering
    └── model2/               Lamar handoff and optional sequence-context extraction
```

Raw FASTQ, BAM, coverage, VCF and large result tables are excluded from version control. The repository contains the analysis logic, frozen summary and provenance needed to regenerate them.
