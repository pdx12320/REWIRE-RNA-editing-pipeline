# REWIRE RNA-editing evidence pipeline

This repository contains the dry-lab workflow used to screen transcriptome-wide C-to-U RNA-editing candidates from three treated and three control RNA-seq libraries.

The pipeline combines:

```text
RNA-seq alignment and preprocessing
→ REDItools2 substitution calling
→ transcript-strand interpretation
→ treated/control comparison
→ exact-allele filtering against a GRCh38-harmonized 293T catalogue
```

![RNA-editing evidence generation pipeline](wiki/assets/figure1_model1_evidence_pipeline.svg)

## Final result

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact 293T catalogue overlaps | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

The retained sites are screening candidates, not a definitive off-target list. The frozen legacy table does not contain independent depth and base counts for control non-calls, and the 293T catalogue is not matched WGS from the exact experimental batch.

## Repository map

| Resource | Purpose |
|---|---|
| [Dry-lab wiki text](wiki/README.md) | Finished workflow, results, contribution and limitations |
| [DBTL record](dbtl/README.md) | Separate development history with four cycles, failed tests and decisions |
| [Pipeline guide](pipeline/README.md) | Reproducible execution steps |
| [Catalogue provenance](pipeline/CATALOGUE_PROVENANCE.md) | Source, liftover and QC |
| [Outputs](pipeline/OUTPUTS.md) | Expected files and frozen result structure |
| [Troubleshooting](pipeline/TROUBLESHOOTING.md) | Observed errors and fixes |
| [Result summary](results/README.md) | Frozen evidence-funnel counts |

## DBTL folder

```text
dbtl/
├── README.md
├── cycle-1-rna-evidence.md
├── cycle-2-public-wgs.md
├── cycle-3-catalogue-harmonization.md
├── cycle-4-final-integration.md
├── failure-log.md
├── decision-log.md
└── reproducibility-checklist.md
```

The wiki intentionally omits the iterative development narrative. The DBTL folder preserves the detailed reasoning and commands without interrupting the public-facing dry-lab story.

## Core implementation

```text
pipeline/scripts/rna/
    RNA processing, substitution calling and evidence integration

pipeline/scripts/catalogue/
    hg18-to-GRCh38 catalogue conversion and exact-allele filtering
```

Large FASTQ, BAM, VCF and per-site result files are excluded from version control. Scripts, QC rules and decision records are retained so that the workflow can be reproduced and audited.
