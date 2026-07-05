# Cycle 5 — Convert screening evidence into LAMAR training labels

## Design

### Question

How can Model 1 provide Model 2 with continuous, background-corrected labels rather than a selected list of caller-positive sites?

### Problem identified

The frozen screening table records whether REDItools2 called a site. Because REDItools2 used an edited-read threshold, a missing control call does not establish zero control editing. The previous strict route measured total depth but did not recover A/C/G/T counts for non-called samples.

## Build

Two scripts were added:

```text
pipeline/scripts/rna/pileup_candidate_base_counts.py
pipeline/scripts/rna/build_lamar_training_table.py
```

The first script directly measures quality-filtered bases at every candidate in all six BAMs. The second script creates transcript-oriented sequence windows and continuous treated-minus-control labels.

A wrapper runs the complete route:

```text
pipeline/scripts/rna/build_lamar_training_labels.sh
```

## Test

The implementation includes unit tests for:

- two-sided Fisher exact calculation;
- Benjamini-Hochberg correction;
- median absolute deviation;
- reverse-complement orientation.

The scripts also stop on missing candidate sites, missing BAMs, ambiguous contig names, invalid sequence-window size or incomplete count tables.

## Learn

### Lesson 1 — Training labels and screening candidates are different products

The 3,333 retained screening candidates are useful for prioritization but should not be the sole Model 2 training set. A predictor trained only on already-filtered positives would inherit the Model 1 selection rules and would have weak negative coverage.

### Lesson 2 — A non-call must be converted into a measured count

Direct pileup allows a control record such as 1 alternate read out of 102 informative reads to remain 0.98% rather than being represented as absent or zero.

### Lesson 3 — Sequence context must follow transcript orientation

Positive-strand genomic C>T and negative-strand genomic G>A events are both converted to a sequence window whose central transcript-oriented base is C. The output includes explicit validation flags so invalid records can be excluded before training.

### Lesson 4 — Statistical significance is not the training target

The principal Model 2 label is the continuous, background-corrected editing efficiency. Pooled Fisher p-values and FDR are retained for screening and QC, not used as substitutes for effect size or replicate-aware validation.

## Next test

Run the route on the six BAM files, inspect label distributions and depth failures, attach PUF-aware metadata, then compare gene-grouped LAMAR fine-tuning against linear, tree-based and one-hot sequence baselines.
