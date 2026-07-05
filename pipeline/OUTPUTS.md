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

Each file should have a `.tbi` index and the expected REDItools2 header.

### Strict one-pass integration outputs

```text
final/CU5.17_EGFP_GC.site_matrix.tsv.gz
final/CU5.17_EGFP_GC.treated_consensus.tsv.gz
final/CU5.17_EGFP_GC.treatment_specific.tsv.gz
```

When the strict route is used, the site matrix retains sample-level call status, independent candidate depth, REDItools2 coverage, alternate-read count, editing fraction and `genomic_catalogue_overlap`.

### Legacy-table compatibility outputs

```text
final_with_293T_catalogue/
├── CU5.17_EGFP_GC.treatment_specific.before_catalogue.annotated.tsv.gz
├── CU5.17_EGFP_GC.treatment_specific.tsv.gz
├── CU5.17_EGFP_GC.catalogue_overlaps.tsv.gz
└── catalogue_integration_summary.tsv
```

The frozen compatibility run contained:

```text
treatment_specific_before_catalogue = 3349
catalogue_overlap                    = 16
final_treatment_specific             = 3333
```

The retained 3,333 rows are catalogue-filtered screening candidates. The 16 excluded rows remain available for audit and must not be treated as negative editing examples.

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

### Integration file used in the frozen run

```text
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
```

### Catalogue QC files

```text
qc/input_header.txt
qc/summary.txt
qc/293T_CG.GRCh38.stats.txt
qc/variants_per_contig.txt
qc/first_20_variants.tsv
qc/checksums.sha256
```

The expected frozen conversion contains 2,885,725 final GRCh38 PASS biallelic SNPs.

## Lamar handoff

```text
$PROJECT/lamar_handoff/
├── CU5.17_EGFP_GC.Lamar_handoff_metadata.tsv
└── CU5.17_EGFP_GC.Lamar_handoff_101nt.tsv   optional sequence-context table
```

The metadata table contains the genomic allele, transcript strand, treated-replicate coverage/alternate counts/edit rates, summary statistics and a provisional `median_treated_edit_rate` label. The sequence-context table additionally contains a fixed-length transcript-oriented sequence with C at the center.
