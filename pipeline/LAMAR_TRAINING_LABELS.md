# Model 1 to Model 2: LAMAR training-label generation

For the production route with complete input manifests, T1 preprocessing audit,
safe checkpoint recovery, direct `samtools mpileup` validation and atomic run
publication, see [Audited Lamar background correction](LAMAR_BACKGROUND_CORRECTION.md).

This module converts the RNA-editing evidence pipeline into a continuous-label dataset for Model 2. It solves the main limitation of a caller-only table: a site that is absent from a REDItools2 output is no longer assumed to have zero alternate reads.

The intended reproducible route preprocesses all samples uniformly. The actual
frozen audit used a Picard MarkDuplicates BAM for T1 and original STAR
coordinate-sorted BAMs for T2/T3/C1/C2/C3. Duplicate-flagged reads were excluded
consistently, but the preprocessing histories were not identical.

## Recommended candidate universe

Use the broadest strand-consistent candidate table available, preferably:

```text
CU5.17_EGFP_GC.site_matrix.tsv.gz
```

Do not train only on the 3,333 final screening candidates. Those records have already passed treated-replicate and control-call filters, so using only that table would create a selected-positive dataset and would not teach LAMAR what low or background editing looks like.

The current 9,930-site matrix is still a called-site universe rather than a complete transcriptome-wide negative set. Model 2 should therefore treat `near_background` records as hard negatives and, when annotation is available, add sequence-matched uncalled editable cytidines as an additional negative set.

The final 3,333 candidates are contained within this broad universe. They cannot
serve as an independent test set for a model trained on the broad matrix.

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

In the frozen audited route, `training_eligible` permits at least 2/3 sufficiently
covered replicates in each group after sequence/orientation QC. High confidence
requires all three treated and all three control replicates plus treated MAD
≤0.05 and control MAD ≤0.02. A missing control measurement never becomes zero;
when group coverage is insufficient, the raw difference and corrected label
remain missing. Corrected labels that genuinely equal zero are retained.

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

Use `prepare_lamar_finetuning_handoff.py` to group overlapping position±50
intervals and exact duplicate sequences before deterministic 80/10/10 assignment.
The optional chromosome strategy is a stronger distribution-shift evaluation.
Do not regenerate a random row-level split from the exported table.

The split file contains all 9,428 eligible rows under one assignment. The
recommended primary subset is the 8,540 high-confidence rows; all eligible rows
and the 7,351 high-confidence, low-control-background rows are sensitivity sets.

Compare LAMAR with at least a linear handcrafted-feature baseline, a tree-based baseline and a one-hot sequence model. Report regression metrics, rank correlation, calibration, performance by editing-rate range and ablation of PUF-aware metadata.

The repository's scalar exporter supervises only the center C. It does not
fabricate a PUF target, binomial total count, or zero-valued labels at other
tokens. `center_index=50` is a nucleotide index before tokenizer special tokens;
validate the tokenizer mapping rather than assuming model index 51. This format
is not claimed to be natively compatible with the historical token-level LAMAR
trainer.

An optional masked token export is available only when the experimentally
confirmed target is supplied explicitly:

```bash
python pipeline/scripts/rna/export_lamar_scalar_regression.py \
  --input /path/to/CU5.17_lamar_splits.tsv.gz \
  --output /path/to/scalar.tsv.gz \
  --token-mask-output /path/to/center_masked.tsv.gz \
  --puf-target-seq "$CONFIRMED_CU517_PUF_TARGET"
```

The mask contains one supervised position at nucleotide index 50. Other label
values are missing/ignored, not known zero, and no total-count vector is emitted.
