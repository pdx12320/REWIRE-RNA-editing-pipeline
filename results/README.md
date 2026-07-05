# Model 1 frozen result summary

The current CU5.17 EGFP-GC analysis produced:

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates with the original depth criteria | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact alleles overlapping the GRCh38-harmonized `293T_CG` catalogue | 16 |
| Final catalogue-filtered candidates | 3,333 |

The catalogue comparison used exact `CHROM:POS:REF:ALT` matching against
`293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz`, which contains 2,885,725 alleles.

Large per-site output tables are not committed here. The compact summary is
stored in `model1_final_summary.tsv`, and the pipeline scripts regenerate the
full tables.
