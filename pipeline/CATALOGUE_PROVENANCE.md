# 293T genomic catalogue provenance

## Source

- Resource: [HEK293 Genome Project v2 data page](https://hek293genome.org/v2/data.php)
- File: `293T_CG.vcf.gz`
- Sample label: 293T
- Variant pipeline: Complete Genomics
- Source assembly: NCBI build 36 / hg18
- Target assembly: GRCh38, using the same FASTA as the RNA-seq branch

The catalogue is an external 293T resource. It is not matched WGS from the CU5.17 experimental cell batch.

## Processing contract

The conversion script performs the following deterministic operations:

1. Detect text or gzip input by file content rather than extension.
2. Retain PASS, biallelic SNPs with a non-reference genotype.
3. Lift coordinates with `hg18ToHg38.over.chain.gz` and CrossMap.
4. Remove records that cannot be mapped.
5. Validate REF alleles against the project GRCh38 FASTA with `bcftools norm -c x`.
6. Normalize, coordinate-sort, bgzip-compress and tabix-index the VCF.
7. Optionally extract a C-to-U-relevant subset containing genomic C→T and G→A alleles.
8. Write summary statistics and SHA-256 checksums.

## Observed QC from the frozen conversion

| Metric | Count |
|---|---:|
| Source PASS biallelic SNPs | 2,914,465 |
| CrossMap-unmapped records | 5,979 |
| GRCh38 REF mismatches removed | 22,761 |
| Final GRCh38 PASS biallelic SNPs | 2,885,725 |
| Final C→T or G→A SNPs | 997,698 |

These counts are an audit target for the specific source file and reference used in this project. A materially different count should trigger inspection of the input VCF, chain file, GRCh38 FASTA and software versions.

## Files available for integration

Full normalized catalogue:

```text
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz.tbi
```

Optional C-to-U-relevant subset:

```text
293T_CG.GRCh38.CtoU_relevant.SNV.vcf.gz
293T_CG.GRCh38.CtoU_relevant.SNV.vcf.gz.tbi
```

The frozen compatibility run used the full normalized catalogue. Because the RNA input table already contained only transcript-oriented C-to-U candidates, exact allele matching against the full catalogue was equivalent for this use case. Sixteen of 3,349 treatment-specific candidates overlapped the catalogue, leaving 3,333 retained screening candidates.

## Interpretation boundary

An exact catalogue match means that the same `CHROM:POS:REF:ALT` allele was reported in the database 293T genome. It supports classification as a plausible genomic variant, but it does not prove that the allele is present in the exact experimental cell batch. Conversely, absence from the catalogue does not prove that a candidate is free of genomic variation.
