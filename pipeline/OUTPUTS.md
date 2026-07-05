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

The site matrix retains sample-level call status, candidate depth, REDItools2 depth, alternate-read count, editing fraction and `genomic_catalogue_overlap`.

## 293T catalogue branch

```text
$CATALOGUE_OUT/
├── 293T_CG.hg18.PASS.biallelic.SNV.vcf.gz
├── 293T_CG.hg18.PASS.biallelic.SNV.vcf.gz.tbi
├── 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
├── 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz.tbi
├── 293T_CG.GRCh38.CtoU_relevant.SNV.vcf.gz
├── 293T_CG.GRCh38.CtoU_relevant.SNV.vcf.gz.tbi
├── logs/
├── qc/
└── tmp/
```

### Recommended integration file

```text
293T_CG.GRCh38.CtoU_relevant.SNV.vcf.gz
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

The expected frozen conversion contains 2,885,725 final GRCh38 PASS biallelic SNPs, including 997,698 C→T or G→A alleles.
