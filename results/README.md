# Legacy ORCA called-site screening summary

This directory preserves the frozen output of the earlier called-site screening analysis:

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact alleles overlapping the GRCh38-harmonized `293T_CG` catalogue | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

The catalogue comparison used exact `CHROM:POS:REF:ALT` matching against `293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz`, which contains 2,885,725 alleles.

The frozen legacy table does not contain independent candidate-site depth and base counts for control non-calls. The 3,333 retained records are therefore an ORCA screening set for experimental prioritization, not a fully depth-qualified off-target set.

These counts are retained for historical audit. They are not the current LAMAR binary positive and negative populations, which are documented in [`LAMAR_BINARY_LABEL_DESIGN.md`](../pipeline/LAMAR_BINARY_LABEL_DESIGN.md).

Large per-site tables are not committed. The compact counts are stored in [`final_summary.tsv`](final_summary.tsv), and the pipeline scripts regenerate the full outputs.

The completed six-sample background-correction audit is summarized in
[`lamar_background_corrected_qc_summary.tsv`](lamar_background_corrected_qc_summary.tsv).
It contains compact broad (9,930-site) and final (3,333-site) QC counts. The full
audit bundle and raw per-sample pileup table remain outside GitHub; only compact
derived model-facing tables are committed below.

The compact derived model handoff in [`../data/processed/`](../data/processed/)
was built without rerunning pileups or changing frozen labels:

| Population | Rows | Intended use |
|---|---:|---|
| All training-eligible | 9,428 | Sensitivity analysis. |
| High confidence | 8,540 | Recommended primary analysis. |
| High confidence, no elevated control flag | 7,351 | Stricter sensitivity analysis; recorded in the full handoff manifest. |
| Corrected-zero among eligible | 1,564 | Retained valid examples. |

The final 3,333 screening candidates are a subset of the broad 9,930 and are not
an independent test set. The broad universe is still candidate-call-derived,
not a complete transcriptome-wide negative universe.
