# Expected outputs

## RNA-seq branch

```text
$PROJECT/
├── fastq/
├── bam/
│   ├── star/
│   ├── markduplicates/
│   └── splitncigarreads/
├── metrics/
├── reditools/
│   ├── coverage/
│   ├── tmp/
│   └── tables/
├── vcf/
├── vep/
├── candidate_depth/          optional strict-depth route
├── final/
├── lamar_training/           Model 1 to Model 2 label route
└── logs/
```

### Per-sample REDItools2 tables

```text
reditools/tables/CU517_GC_T1.txt.gz
reditools/tables/CU517_GC_T2.txt.gz
reditools/tables/CU517_GC_T3.txt.gz
reditools/tables/CU517_GC_C1.txt.gz
reditools/tables/CU517_GC_C2.txt.gz
reditools/tables/CU517_GC_C3.txt.gz
```

Each file should have the expected REDItools2 header and a tabix index when required by downstream tools.

## Strict integration outputs

```text
final/CU5.17_EGFP_GC.site_matrix.tsv.gz
final/CU5.17_EGFP_GC.treated_consensus.tsv.gz
final/CU5.17_EGFP_GC.treatment_specific.tsv.gz
```

The strict route retains sample-level call status, independent candidate depth, REDItools2 coverage, alternate-read count, editing fraction and `genomic_catalogue_overlap`.

## LAMAR training-label outputs

```text
lamar_training/
├── base_counts/
│   ├── CU517_GC_T1.candidate_base_counts.tsv.gz
│   ├── CU517_GC_T2.candidate_base_counts.tsv.gz
│   ├── CU517_GC_T3.candidate_base_counts.tsv.gz
│   ├── CU517_GC_C1.candidate_base_counts.tsv.gz
│   ├── CU517_GC_C2.candidate_base_counts.tsv.gz
│   └── CU517_GC_C3.candidate_base_counts.tsv.gz
└── CU5.17_EGFP_GC.lamar_training_labels.tsv.gz
```

Each per-sample table contains quality-filtered A/C/G/T counts, forward/reverse counts, ref/alt counts, allele depth, editing rate and exclusion diagnostics. The combined table retains all six sample measurements and adds transcript-oriented sequence windows, treated/control medians, corrected editing efficiency, replicate MAD, pooled Fisher/FDR, label class, confidence and training eligibility.

Prefer the broad strand-consistent site matrix as input. The final 3,333-site screening table is already selected for treated presence and control non-calls and should not be the sole Model 2 training set.

## Legacy-table compatibility outputs

```text
final_with_293T_catalogue/
├── CU5.17_EGFP_GC.treatment_specific.before_catalogue.annotated.tsv.gz
├── CU5.17_EGFP_GC.treatment_specific.tsv.gz
├── CU5.17_EGFP_GC.catalogue_overlaps.tsv.gz
└── catalogue_integration_summary.tsv
```

Frozen result:

```text
treatment-specific before catalogue    3,349
exact catalogue overlaps                  16
retained screening candidates          3,333
```

The 16 excluded alleles are preserved for audit. They should not be treated as negative editing examples.

## 293T catalogue branch

```text
$CATALOGUE_OUT/
├── 293T_CG.hg18.PASS.biallelic.SNV.vcf.gz
├── 293T_CG.hg18.PASS.biallelic.SNV.vcf.gz.tbi
├── 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
├── 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz.tbi
├── logs/
├── qc/
└── tmp/
```

The frozen integration used:

```text
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
```

Catalogue QC files:

```text
qc/input_header.txt
qc/summary.txt
qc/293T_CG.GRCh38.stats.txt
qc/variants_per_contig.txt
qc/first_20_variants.tsv
qc/checksums.sha256
```

The frozen conversion contains 2,885,725 final GRCh38 PASS biallelic SNPs.
