# Data policy and compact model handoff

`data/processed/` contains only compact, derived, model-facing products from the
completed CU5.17 six-sample audit. It does not contain BAM, FASTQ, reference
FASTA, full pileup tables, absolute symlinks or runtime configuration tied to a
user-specific server path.

| File | Purpose |
|---|---|
| `CU5.17_lamar_all_eligible.tsv.gz` | All 9,428 frozen training-eligible rows; sensitivity analysis. |
| `CU5.17_lamar_high_confidence.tsv.gz` | 8,540-row recommended primary dataset. |
| `CU5.17_lamar_splits.tsv.gz` | Deterministic 80/10/10 assignment for all eligible rows with overlap and exact-sequence group IDs. |
| `data_dictionary.tsv` | Column types and definitions. |
| `handoff_manifest.json` | Input hashes, frozen formula, counts, split algorithm and guardrails. |
| `checksums.sha256` | SHA-256 checksums for the committed compact files. |

The complete handoff builder also produces the stricter high-confidence,
low-control-background sensitivity table, the excluded table, `split_qc.json`
and a generated README in the requested output directory.

The corrected target is
`max(median(treated editing rates) - median(control editing rates), 0)`. Missing
control evidence is never filled with zero, while valid corrected-zero examples
are retained. The final 3,333 candidates are a subset of the broad 9,930-site
matrix and are not an independent model test set.

The broad matrix is derived from called candidate sites and is not a complete
transcriptome-wide negative universe. Sequence-matched, sufficiently covered
uncalled cytidines remain a future improvement. These derived outputs passed
computational QC; they are not experimentally or biologically validated labels.

Regenerate the public files from the frozen audited tables without touching BAMs:

```bash
python pipeline/scripts/rna/prepare_lamar_finetuning_handoff.py \
  --labels /path/to/background_corrected_labels.tsv.gz \
  --metadata /path/to/lamar_ready_metadata.tsv.gz \
  --output-dir /path/to/full_handoff \
  --public-copy-dir data/processed \
  --seed 20260715 \
  --split-strategy overlap_cluster
```

Do not add `puf_target_seq` until the wet-lab/model team supplies the exact
CU5.17 target. Do not infer `label_total_count` from the median-difference label.
