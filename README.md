# REWIRE RNA-editing evidence pipeline

ORCA is a programmable PUF–APOBEC RNA-editing system. This repository documents how six RNA-seq libraries were converted into auditable sequence-level labels for LAMAR candidate prioritization.

The current binary dataset does not define negatives as sites omitted by an editing caller. It directly measures expressed transcript cytosines and requires complete six-sample evidence before assigning a strict computational-negative label.

![Current computational-positive and strict computational-negative label design](wiki/assets/figure1_current_binary_label_design.svg)

## Current model-facing dataset

| Label population | Sites | Interpretation |
|---|---:|---|
| Computational positives | **1,513** | Corrected editing efficiency greater than 0.10 after coverage, control-background, sequence, complexity, and WGS checks |
| High-confidence positive audit subset | **1,457** | Main positives with complete six-sample coverage and stronger replicate-consistency criteria |
| Strict computational negatives | **2,821,734** | Expressed exon C sites with depth at least 20 and zero target-ALT reads in all six samples |

Every model sequence is transcript-oriented, contains 101 nucleotides, and has C at zero-based index 50. Coverage and annotation were used to construct and audit labels, not as default model inputs.

## Evidence design

- The same pileup implementation was used for positive and negative recounting.
- Computational positives required coverage in at least two treated and two control replicates, corrected efficiency above 0.10, and control median at most 0.02.
- Strict computational negatives required usable depth of at least 20 and target-ALT count equal to zero in every treated and control replicate.
- Both classes excluded central WGS variants, invalid sequence orientation, incomplete windows, ambiguous centers, and low-complexity contexts.
- All original broad candidates and all positive definitions were excluded from the strict negative universe.
- Gene, genomic-center, overlapping-window, and exact-sequence relations were grouped before splitting.

The detailed frozen specification is available in [LAMAR binary label design](pipeline/LAMAR_BINARY_LABEL_DESIGN.md).

## Legacy screening analysis

The earlier called-site screening funnel is no longer presented as the model-facing dataset. It did not construct a transcriptome-wide, fully depth-qualified negative class.

The legacy screening counts remain in [results](results/README.md) and the [DBTL record](dbtl/README.md) for historical audit and reproducibility. They must not be interpreted as the current positive-to-negative label pipeline.

## Model interface

The current binary route measures quality-filtered A, C, G, and T counts in all six MarkDuplicates BAMs. It preserves replicate-level depth and target-ALT evidence, then extracts a 101-nucleotide transcript-oriented sequence.

See [LAMAR training-label generation](pipeline/LAMAR_TRAINING_LABELS.md) for the legacy continuous-label route and [audited background correction](pipeline/LAMAR_BACKGROUND_CORRECTION.md) for its frozen QC record.

## Repository map

| Resource | Purpose |
|---|---|
| [Dry-lab wiki](wiki/README.md) | Finished workflow, results, contribution and limitations |
| [DBTL record](dbtl/README.md) | Development cycles plus failure and decision logs |
| [Pipeline guide](pipeline/README.md) | Reproducible commands |
| [LAMAR binary label design](pipeline/LAMAR_BINARY_LABEL_DESIGN.md) | Current computational-positive and strict computational-negative specification |
| [Legacy LAMAR label guide](pipeline/LAMAR_TRAINING_LABELS.md) | Earlier continuous-label construction route |
| [Catalogue provenance](pipeline/CATALOGUE_PROVENANCE.md) | Source, liftover and quality control |
| [Outputs](pipeline/OUTPUTS.md) | Expected files |
| [Troubleshooting](pipeline/TROUBLESHOOTING.md) | Observed errors and fixes |
| [Legacy result summary](results/README.md) | Frozen historical screening counts |

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
