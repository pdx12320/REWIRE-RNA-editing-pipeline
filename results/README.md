# ORCA frozen result summary

The current CU5.17 EGFP-GC analysis produced:

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact alleles overlapping the GRCh38-harmonized `293T_CG` catalogue | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

The catalogue comparison used exact `CHROM:POS:REF:ALT` matching against `293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz`, which contains 2,885,725 alleles.

The frozen legacy table does not contain independent candidate-site depth and base counts for control non-calls. The 3,333 retained records are therefore an ORCA screening set for experimental prioritization, not a fully depth-qualified off-target set.

Large per-site tables are not committed. The compact counts are stored in [`final_summary.tsv`](final_summary.tsv), and the pipeline scripts regenerate the full outputs.

The completed six-sample background-correction audit is summarized in
[`lamar_background_corrected_qc_summary.tsv`](lamar_background_corrected_qc_summary.tsv).
It contains compact broad (9,930-site) and final (3,333-site) QC counts; the large
per-site labels and pileups remain outside GitHub.
