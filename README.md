# REWIRE RNA-editing evidence pipeline

This repository contains the dry-lab workflow used to screen transcriptome-wide C-to-U RNA-editing candidates from three treated and three control RNA-seq libraries.

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

## Repository map

| Resource | Purpose |
|---|---|
| [Dry-lab wiki](wiki/README.md) | Finished workflow, results, contribution and limitations |
| [DBTL record](dbtl/README.md) | Ten development cycles plus failure and decision logs |
| [Pipeline guide](pipeline/README.md) | Reproducible commands |
| [Catalogue provenance](pipeline/CATALOGUE_PROVENANCE.md) | Source, liftover and QC |
| [Outputs](pipeline/OUTPUTS.md) | Expected files |
| [Troubleshooting](pipeline/TROUBLESHOOTING.md) | Observed errors and fixes |
| [Result summary](results/README.md) | Frozen counts |

## DBTL cycles

```text
1  RNA evidence
2  public WGS test
3  catalogue harmonization
4  final integration
5  GATK read groups
6  REDItools2 environment
7  contigs and coverage
8  control-depth evidence
9  legacy compatibility
10 audit and reporting
```

The retained sites are screening candidates. The frozen legacy table does not contain independent depth and base counts for control non-calls, and the 293T catalogue is external to the exact experimental batch.

Large sequencing and intermediate files remain outside GitHub. The repository retains scripts, environment definitions, QC rules, frozen counts and decision records.
