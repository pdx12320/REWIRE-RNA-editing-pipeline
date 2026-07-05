# Cycle 3 — Harmonize the HEK293 293T_CG catalogue

## Design

### Question

Could the database-released `293T_CG` variant catalogue provide a more credible genomic exclusion resource than the failed public-WGS reconstruction?

### Source

The selected resource was the HEK293 Genome Project `293T_CG.vcf.gz`, generated with the Complete Genomics pipeline.

Important constraints:

```text
source assembly: NCBI build 36 / hg18
RNA-seq assembly: GRCh38
sample relation: external 293T catalogue, not matched WGS
```

The catalogue could only be compared with RNA candidates after assembly conversion and REF-allele validation.

### Design requirements

The conversion had to:

1. detect whether the input was text or gzip-compressed;
2. retain PASS biallelic SNPs with a non-reference genotype;
3. convert coordinates from hg18 to GRCh38;
4. remove unmapped records;
5. verify each REF allele against the same GRCh38 FASTA used by the RNA branch;
6. normalize and sort the VCF;
7. bgzip-compress and tabix-index the output;
8. record QC counts and checksums.

## Build

### Environment

```bash
conda env create -f pipeline/env/genomic_catalogue.yml
conda activate rewire_catalogue
```

The environment provides CrossMap, bcftools, samtools, bgzip and tabix.

### Inputs

```bash
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
CATALOGUE_IN=/data/ydx/igem/293T_CG.vcf
CHAIN=/data/ydx/igem/hg18ToHg38.over.chain.gz
CATALOGUE_OUT=/data/ydx/igem/293T_CG_GRCh38
```

### Conversion command

```bash
bash pipeline/scripts/catalogue/process_293T_CG_to_GRCh38.sh \
  --input "$CATALOGUE_IN" \
  --reference "$REF" \
  --chain "$CHAIN" \
  --outdir "$CATALOGUE_OUT" \
  --threads 16
```

### Implemented data flow

```text
input-format normalization
→ PASS biallelic SNP selection
→ CrossMap hg18-to-GRCh38 liftover
→ remove unmapped records
→ GRCh38 REF validation
→ remove REF mismatches
→ normalize and coordinate-sort
→ bgzip compression
→ tabix indexing
→ QC summary and SHA-256 checksums
```

### Principal output

```text
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz.tbi
```

An optional C-to-U-relevant C>T/G>A subset can also be generated, but the frozen final comparison used the full normalized catalogue because the RNA input was already restricted to transcript-oriented C-to-U candidates.

## Test

### Frozen QC counts

| Stage | Count |
|---|---:|
| Source PASS biallelic SNPs | 2,914,465 |
| CrossMap-unmapped records | 5,979 |
| GRCh38 REF mismatches removed | 22,761 |
| Final GRCh38 PASS biallelic SNPs | 2,885,725 |
| C>T or G>A subset | 997,698 |

The final full catalogue contained 2,885,725 exact alleles and was readable by bcftools and tabix.

### Validation checks

```bash
gzip -t 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
bcftools index -n 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
bcftools stats 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
```

REF validation used the same GRCh38 FASTA as the RNA workflow. This prevented comparison of lifted coordinates carrying an incompatible reference allele.

## Learn

### Failure 1 — Misleading file extension

Observed symptom:

```text
Exec format error
```

Diagnosis: the source file content was gzip-compressed even when the local filename ended in `.vcf`.

Change: the script now detects compression from file bytes rather than trusting the extension.

### Failure 2 — Chain download failure

Observed symptom: the UCSC chain download timed out on the server.

Change: the workflow accepts a validated local copy of `hg18ToHg38.over.chain.gz` and checks it with `gzip -t`.

### Failure 3 — Unsorted CrossMap output

Observed symptom:

```text
Unsorted positions on sequence ...
index: failed to create index
```

Diagnosis: coordinate liftover did not guarantee final sort order.

Change: `bcftools sort` is run before bgzip/tabix indexing.

### Failure 4 — REF mismatches after liftover

Observed result: 22,761 lifted records did not match the target GRCh38 reference allele.

Diagnosis: a mapped coordinate alone is insufficient; the allele must also be compatible with the target reference.

Change:

```bash
bcftools norm -f GRCh38.fa -c x
```

is used to remove incompatible records.

### Final decision

The harmonized `293T_CG` catalogue replaced the failed public-WGS reconstruction because it provided a genome-scale, QC-audited external variant resource.

It remained an exclusion catalogue rather than definitive matched-genome evidence because it was generated from another 293T sample and calling pipeline.

## Output of this cycle

```text
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz.tbi
logs/
qc/
tmp/
```

The full catalogue became the genomic input for exact-allele integration in Cycle 4.
