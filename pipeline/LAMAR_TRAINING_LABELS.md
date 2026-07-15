# Model 1 to Model 2: LAMAR training-label generation

For the production route with complete input manifests, T1 preprocessing audit,
safe checkpoint recovery, direct `samtools mpileup` validation and atomic run
publication, see [Audited Lamar background correction](LAMAR_BACKGROUND_CORRECTION.md).

This module converts the RNA-editing evidence pipeline into a continuous-label dataset for Model 2. It solves the main limitation of a caller-only table: a site that is absent from a REDItools2 output is no longer assumed to have zero alternate reads.

## Recommended candidate universe

Use the broadest strand-consistent candidate table available, preferably:

```text
CU5.17_EGFP_GC.site_matrix.tsv.gz
```

Do not train only on the 3,333 final screening candidates. Those records have already passed treated-replicate and control-call filters, so using only that table would create a selected-positive dataset and would not teach LAMAR what low or background editing looks like.

The current 9,930-site matrix is still a called-site universe rather than a complete transcriptome-wide negative set. Model 2 should therefore treat `near_background` records as hard negatives and, when annotation is available, add sequence-matched uncalled editable cytidines as an additional negative set.

## Stage 1: quality-filtered base counts in all six BAMs

`pileup_candidate_base_counts.py` measures every candidate directly in every treated and control BAM. It records:

- A, C, G and T counts;
- forward and reverse counts for each base;
- reference and alternate counts;
- allele depth and editing rate;
- alternate-strand balance;
- observations excluded for low mapping quality, low base quality, duplicate, secondary, supplementary, QC-fail, deletion/ref-skip or non-ACGT bases.

Default filters are mapping quality 30 and base quality 20. Counts are produced even when the alternate allele is below the REDItools2 `-me 20` calling threshold.

## Stage 2: continuous LAMAR labels

`build_lamar_training_table.py` combines the six count tables and creates a 101-nt transcript-oriented sequence window centered on the editable cytidine.

For replicate editing rates `r`, the principal label is:

```text
corrected_editing_efficiency = max(
    0,
    median(treated editing rates) - median(control editing rates)
)
```

The table also retains:

- treated and control median editing rates;
- pooled editing rates;
- replicate median absolute deviations;
- raw treated-minus-control difference;
- pooled Fisher exact p-value and BH-FDR;
- depth coverage in each group;
- sequence and strand validation flags;
- `training_eligible`, `label_confidence` and `label_class`.

The pooled Fisher test is a screening statistic based on read counts. It is not a replacement for replicate-aware biological inference, because reads are not independent biological replicates.

## Run the complete label-generation route

```bash
conda env create -f pipeline/env/lamar_labels.yml
conda activate rewire_lamar_labels

PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
CANDIDATES="$PROJECT/final/CU5.17_EGFP_GC.site_matrix.tsv.gz"

bash pipeline/scripts/rna/build_lamar_training_labels.sh \
  "$PROJECT" \
  "$REF" \
  "$CANDIDATES" \
  pipeline/config/samples.tsv \
  pipeline/config/site_metadata.tsv
```

The optional site-metadata table is joined by exact `CHROM:POS:REF:ALT`. It can provide `gene_id`, `transcript_id`, `puf_target_sequence`, PUF similarity, editor-to-target distance and region annotations. An example is provided in `pipeline/config/site_metadata.example.tsv`.

## Output

```text
$PROJECT/lamar_training/
├── base_counts/
│   ├── CU517_GC_T1.candidate_base_counts.tsv.gz
│   ├── CU517_GC_T2.candidate_base_counts.tsv.gz
│   ├── CU517_GC_T3.candidate_base_counts.tsv.gz
│   ├── CU517_GC_C1.candidate_base_counts.tsv.gz
│   ├── CU517_GC_C2.candidate_base_counts.tsv.gz
│   └── CU517_GC_C3.candidate_base_counts.tsv.gz
└── CU5.17_EGFP_GC.lamar_training_labels.tsv.gz
```

## Model 2 split and evaluation rules

Model 2 should not use a random site-level split as its primary result. Nearby sites and sites from the same transcript can share highly similar sequence windows.

Use, in order of preference:

1. gene-grouped or transcript-grouped cross-validation;
2. held-out PUF target or experiment when multiple targets become available;
3. independent experimental validation.

Compare LAMAR with at least a linear handcrafted-feature baseline, a tree-based baseline and a one-hot sequence model. Report regression metrics, rank correlation, calibration, performance by editing-rate range and ablation of PUF-aware metadata.
