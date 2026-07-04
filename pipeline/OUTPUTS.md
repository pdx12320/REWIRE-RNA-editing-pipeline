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
├── candidate_depth/
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

Each file should have a `.tbi` index and the 14-column REDItools2 header.

### Final RNA evidence tables

```text
final/CU5.17_EGFP_GC.site_matrix.tsv.gz
final/CU5.17_EGFP_GC.treated_consensus.tsv.gz
final/CU5.17_EGFP_GC.treatment_specific.tsv.gz
```

The site matrix retains sample-level call status, candidate depth, REDItools2 depth, alternate-read count, editing fraction and WGS overlap.

## WGS branch

```text
$WGS_OUT/
├── metadata/
├── sra/
├── fastq/
├── bam/
├── qc/
├── vcf/
├── logs/
└── tmp/
```

### Merge mode

```text
vcf/HEK293T_3runs.filtered.SNV.vcf.gz
```

### Consensus mode

```text
vcf/HEK293T_3runs.consensus2of3.SNV.vcf.gz
vcf/HEK293T_3runs.union.SNV.vcf.gz
```

The two-of-three consensus is the conservative exclusion blacklist. The union is better treated as a broad annotation flag.
