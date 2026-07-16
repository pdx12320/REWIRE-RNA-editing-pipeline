# REWIRE RNA-editing pipeline

[![LAMAR label tests](https://github.com/pdx12320/REWIRE-RNA-editing-pipeline/actions/workflows/lamar-label-tests.yml/badge.svg)](https://github.com/pdx12320/REWIRE-RNA-editing-pipeline/actions/workflows/lamar-label-tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11–3.12](https://img.shields.io/badge/Python-3.11%20%7C%203.12-blue.svg)](pipeline/env/lamar_labels.yml)

ORCA is the programmable PUF–APOBEC RNA-editing system evaluated by this repository. The dry-lab pipeline screens transcriptome-wide C-to-U candidates from three treated and three control RNA-seq libraries. The repository also contains the Model 1 interface that converts this evidence into continuous labels for downstream LAMAR training.

The repository has three connected but distinct layers:

| Layer | Purpose |
|---|---|
| Transcriptome-wide screening | Align RNA-seq, call substitutions, interpret transcript strand and compare exact alleles with a harmonized 293T catalogue. |
| Background-corrected labels | Recount every candidate in six BAMs and compute audited replicate-median treated-minus-control targets. |
| LAMAR fine-tuning handoff | Validate the frozen tables, build leakage-resistant splits and export center-C scalar-regression records. |

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

The frozen audited run did **not** use identical preprocessing histories. T1 used
the available coordinate-sorted Picard MarkDuplicates BAM; T2, T3, C1, C2 and C3
used original STAR coordinate-sorted BAMs. Duplicate-flagged reads were excluded
with the same pileup flag filter in all samples, but the five STAR BAMs had never
been duplicate-marked. This exception remains part of the QC record.

See [LAMAR training-label generation](pipeline/LAMAR_TRAINING_LABELS.md) for the
composable route and [audited background correction](pipeline/LAMAR_BACKGROUND_CORRECTION.md)
for the production six-sample run, recovery checks, frozen QC and scalar-versus-token
fine-tuning boundary.

## Model teammate quick start

The handoff builder is stdlib-only. The optional baseline environment adds
scikit-learn and uses Python 3.12; CI also tests the repository on Python 3.11.

```bash
conda env create -f pipeline/env/lamar_scalar_baseline.yml
conda activate rewire_lamar_scalar

python -m unittest discover -s tests -p "test_*.py" -v

# Inspect the committed production split QC embedded in the manifest.
python -m json.tool data/processed/handoff_manifest.json

# Run directly on the committed recommended-primary split assignments.
python pipeline/scripts/rna/export_lamar_scalar_regression.py \
  --input data/processed/CU5.17_lamar_splits.tsv.gz \
  --output /tmp/CU5.17_lamar_scalar_high_confidence.tsv.gz

python examples/train_scalar_baseline.py \
  --input data/processed/CU5.17_lamar_splits.tsv.gz \
  --subset high_confidence \
  --output-json /tmp/CU5.17_scalar_baseline_metrics.json

# Rebuild a full handoff later from the two frozen audit inputs.

LABELS=/path/to/background_corrected_labels.tsv.gz
METADATA=/path/to/lamar_ready_metadata.tsv.gz
HANDOFF=/path/to/CU5.17_lamar_finetuning_handoff

python pipeline/scripts/rna/prepare_lamar_finetuning_handoff.py \
  --labels "$LABELS" \
  --metadata "$METADATA" \
  --output-dir "$HANDOFF" \
  --seed 20260715 \
  --split-strategy overlap_cluster

python -m json.tool "$HANDOFF/split_qc.json"

python pipeline/scripts/rna/export_lamar_scalar_regression.py \
  --input "$HANDOFF/CU5.17_lamar_splits.tsv.gz" \
  --output "$HANDOFF/CU5.17_lamar_scalar_high_confidence.tsv.gz"

python examples/train_scalar_baseline.py \
  --input "$HANDOFF/CU5.17_lamar_splits.tsv.gz" \
  --subset high_confidence \
  --output-json "$HANDOFF/scalar_baseline_metrics.json"
```

Do not add a random row-level split after export. `center_index=50` is the
zero-based nucleotide position before tokenization; verify any shift introduced
by `[CLS]` or other tokenizer special tokens instead of hard-coding model token
index 51.

## Model-facing datasets

| Dataset or population | Rows | Recommended use |
|---|---:|---|
| Broad called-candidate universe | 9,930 | Audited source universe; not a complete transcriptome-wide negative universe. |
| All training-eligible | 9,428 | Sensitivity analysis; eligibility allows at least 2/3 sufficiently covered replicates per group. |
| High confidence | 8,540 | **Recommended primary fine-tuning dataset**; requires all 3+3 covered replicates and frozen replicate-consistency thresholds. |
| High confidence, low control background | 7,351 | Stricter sensitivity analysis. |
| Zero corrected labels among eligible rows | 1,564 | Keep as valid background-corrected examples; do not delete or relabel. |
| Final selected candidates | 3,333 | Screening/prioritization subset of the broad matrix; never use as an independent test set against broad-matrix training. |

The broad 9,930-site universe is derived from called candidate sites. A future
negative-set improvement is to add sequence-matched, sufficiently covered,
uncalled transcript-oriented cytidines. The current handoff does not invent
those sites, `puf_target_seq`, `label_total_count`, or non-center token labels.

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
| [Compact model data](data/README.md) | Public derived tables, manifests and checksums |

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
