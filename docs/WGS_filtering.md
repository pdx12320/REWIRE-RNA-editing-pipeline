# Public HEK293T WGS processing and genomic-SNV filtering

## Purpose

RNA-seq mismatches are not automatically RNA-editing events. A C-to-T or G-to-A signal may also reflect a genomic single-nucleotide variant already present in the HEK293T genome. This module builds an external HEK293T genomic-variant blacklist from three public WGS runs and flags RNA candidates that match the same `CHROM:POS:REF:ALT` allele.

The public runs configured in this repository are:

```text
SRR37832939
SRR37832940
SRR37832941
```

Because these runs were not generated from the exact CU5.17 experimental cell batch, the output is described as an **external public HEK293T blacklist**, not matched WGS.

## Software

- SRA Toolkit: download and FASTQ conversion
- BWA-MEM2 or BWA-MEM: short-read alignment to GRCh38
- samtools: sorting, indexing and BAM validation
- GATK MarkDuplicates: duplicate marking
- bcftools mpileup/call: SNV discovery
- custom Python scripts: per-run filtering and N-of-M consensus construction

## 1. Create the environment

```bash
conda env create -f environment/wgs_pipeline.yml
conda activate rewire_wgs
```

Check:

```bash
prefetch --version
bwa-mem2 version 2>/dev/null || bwa 2>&1 | head
samtools --version
bcftools --version
gatk --version
```

## 2. Inspect SRA metadata before processing

Consecutive SRR accessions do not prove that runs represent one biological sample. The metadata script checks library strategy, layout, organism and BioSample accession:

```bash
bash scripts/wgs/00_check_sra_metadata.sh \
  config/wgs_runs.tsv \
  wgs_sra_metadata.tsv
```

Interpretation:

- same `sample_accession`: use merge mode;
- different `sample_accession`: call each run separately and retain exact alleles supported by at least two runs;
- all runs must be paired-end WGS.

## 3. Run the WGS pipeline

```bash
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
OUT=/data/ydx/igem/HEK293T_public_WGS_3runs

nohup bash scripts/wgs/run_3run_wgs_pipeline.sh \
  --runs config/wgs_runs.tsv \
  --reference "$REF" \
  --outdir "$OUT" \
  --mode auto \
  --threads 32 \
  --min-dp 10 \
  --min-alt 3 \
  --min-vaf 0.05 \
  --min-qual 20 \
  > "$OUT.pipeline.log" 2>&1 &
```

The pipeline performs:

```text
SRA metadata validation
→ prefetch and fasterq-dump
→ BWA-MEM2/BWA-MEM alignment to the same GRCh38 reference used for RNA-seq
→ coordinate sorting and BAM indexing
→ GATK MarkDuplicates
→ bcftools mpileup and multiallelic calling
→ reference normalization and SNP-only filtering
→ merge-mode VCF or 2-of-3 exact-allele consensus blacklist
```

## 4. Default WGS thresholds

A variant is retained within a single WGS call set when:

```text
site depth ≥10
alternate reads ≥3
alternate-allele fraction ≥0.05
variant QUAL ≥20
FILTER is PASS or unset
```

These thresholds are intended to flag plausible genomic alleles rather than produce a clinical germline call set.

## 5. Outputs

When all three runs belong to one BioSample:

```text
$OUT/vcf/HEK293T_3runs.filtered.SNV.vcf.gz
```

When runs represent distinct BioSamples:

```text
$OUT/vcf/HEK293T_3runs.consensus2of3.SNV.vcf.gz
$OUT/vcf/HEK293T_3runs.union.SNV.vcf.gz
```

Use the 2-of-3 consensus as the conservative exclusion blacklist. The union is broader and is better used as an annotation flag.

## 6. Apply the blacklist to Model 1

```bash
python3 scripts/filter_c_to_u_and_compare.py \
  --manifest config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --wgs-vcf "$OUT/vcf/HEK293T_3runs.consensus2of3.SNV.vcf.gz" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

The existing filter compares exact chromosome, position, reference and alternate alleles. A coordinate overlap with a different alternate base is not removed.

## Interpretation boundary

This public blacklist strengthens genomic-SNV filtering but does not fully substitute for WGS from the exact experimental HEK293T batch. The final output should therefore be described as **treatment-associated RNA-editing candidates filtered against an external HEK293T genomic-variant catalogue**.
