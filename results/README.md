# Model 1 frozen result summary

The current CU5.17 EGFP-GC analysis produced:

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact alleles overlapping the GRCh38-harmonized `293T_CG` catalogue | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

The catalogue comparison used exact `CHROM:POS:REF:ALT` matching against
`293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz`, which contains 2,885,725 alleles.

The frozen legacy result table does not contain independent candidate-site
depth and base counts for control non-calls. The 3,333 retained records are
therefore a screening set for Lamar inference and experimental prioritization,
not a fully depth-qualified off-target set.

Large per-site output tables are not committed here. The compact summary is
stored in `model1_final_summary.tsv`, and the pipeline scripts regenerate the
full tables.
