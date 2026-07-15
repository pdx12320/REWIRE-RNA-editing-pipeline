# Transcriptome-wide RNA-editing evidence pipeline

ORCA is the programmable PUF–APOBEC RNA-editing system evaluated by this repository. The dry-lab pipeline screens transcriptome-wide C-to-U candidates from three treated and three control RNA-seq libraries. The repository also contains the Model 1 interface that converts this evidence into continuous labels for downstream LAMAR training.

```text
RNA-seq alignment and preprocessing
→ REDItools2 substitution calling
→ transcript-strand interpretation
→ treated/control comparison
→ exact-allele filtering against a GRCh38-harmonized 293T catalogue
→ six-sample candidate base counting
→ background-corrected LAMAR training labels
```

![Transcriptome-wide RNA-editing screening workflow](wiki/assets/figure1_screening_workflow.svg)

## Final result

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact 293T catalogue overlaps | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

The retained sites are screening candidates rather than confirmed off-targets. The frozen legacy table does not contain independent depth and base counts for control non-calls, and the 293T catalogue is external to the exact experimental batch.

## Model 1 to Model 2 interface

The LAMAR label route measures quality-filtered A/C/G/T counts at every candidate in all six BAMs, preserves replicate-level ref/alt counts, extracts transcript-oriented sequence windows and calculates background-corrected editing efficiency. The recommended training universe is the broad 9,930-site matrix rather than only the 3,333 already-filtered candidates.

See [LAMAR training-label generation](pipeline/LAMAR_TRAINING_LABELS.md) for the
composable route and [audited background correction](pipeline/LAMAR_BACKGROUND_CORRECTION.md)
for the production six-sample run, recovery checks, frozen QC and scalar-versus-token
fine-tuning boundary.

## Repository map

| Resource | Purpose |
|---|---|
| [Dry-lab wiki](wiki/README.md) | Finished workflow, results, contribution and limitations |
| [DBTL record](dbtl/README.md) | Development cycles plus failure and decision logs |
| [Pipeline guide](pipeline/README.md) | Reproducible commands |
| [LAMAR label guide](pipeline/LAMAR_TRAINING_LABELS.md) | Six-sample base counting and Model 2 training-table construction |
| [Catalogue provenance](pipeline/CATALOGUE_PROVENANCE.md) | Source, liftover and quality control |
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
11 LAMAR training labels
12 audited six-sample background correction
```

Large sequencing and intermediate files remain outside GitHub. The repository retains scripts, environment definitions, quality-control rules, frozen counts and decision records.
