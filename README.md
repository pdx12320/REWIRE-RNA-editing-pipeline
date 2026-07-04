# REWIRE Model 1 — RNA-editing evidence pipeline

This repository contains the complete implementation of **Model 1**, which converts treated/control RNA-seq and public HEK293T WGS data into an auditable evidence matrix for candidate transcript-level C-to-U editing.

## Main links

- **Main iGEM Wiki page:** [`wiki/Model1_RNA_editing_evidence_pipeline.md`](wiki/Model1_RNA_editing_evidence_pipeline.md)
- **Copy-ready Methods section:** [`wiki/Model1_Methods_copy_ready.md`](wiki/Model1_Methods_copy_ready.md)
- **Figure 1 with RNA-seq and WGS branches:** [`wiki/assets/model1_pipeline_with_wgs.svg`](wiki/assets/model1_pipeline_with_wgs.svg)
- **RNA sample manifest:** [`config/samples.tsv`](config/samples.tsv)
- **Public WGS run manifest:** [`config/wgs_runs.tsv`](config/wgs_runs.tsv)
- **WGS technical guide:** [`docs/WGS_filtering.md`](docs/WGS_filtering.md)
- **Executable scripts:** [`scripts/`](scripts/)

## Scientific question

Reporter editing demonstrates on-target activity but does not establish transcriptome-wide specificity. Model 1 asks which substitutions are:

1. reproducible across three treated RNA-seq replicates;
2. absent or minimal in three controls with sufficient coverage;
3. consistent with transcript-level C-to-U editing after strand annotation;
4. not readily explained by a public HEK293T genomic SNV.

## Figure 1

![RNA-editing evidence generation pipeline](wiki/assets/model1_pipeline_with_wgs.svg)

The RNA-seq and WGS branches use the same GRCh38 reference and converge in an exact-allele, site-level evidence matrix. Public WGS is an external blacklist rather than WGS matched to the experimental CU5.17 cell batch.

## Workflow

```text
RNA-seq branch
SRA → FASTQ → STAR → GATK → coverage map → REDItools2 → VEP
      → all-sample depth → treated/control comparison

WGS branch
SRA metadata → FASTQ → BWA-MEM2/BWA-MEM → GATK MarkDuplicates
      → bcftools calling → merged or 2-of-3 consensus SNV blacklist

Integration
RNA evidence matrix + exact-allele WGS blacklist
      → treatment-associated C-to-U candidates
```

## Repository layout

```text
config/
  samples.tsv
  wgs_runs.tsv

environment/
  reditools2_py2.yml
  wgs_pipeline.yml

scripts/
  download_sra_fastq.py
  run_star_alignment.py
  run_gatk_preprocessing.py
  generate_reditools_coverage_limited.sh
  run_reditools_all_samples.sh
  reditools_union_to_vcf.py
  run_vep_annotation.py
  build_candidate_depth_tables.sh
  filter_c_to_u_and_compare.py
  wgs/
    00_check_sra_metadata.sh
    run_3run_wgs_pipeline.sh
    filter_single_sample_vcf.py
    build_consensus_blacklist.py

wiki/
  Model1_RNA_editing_evidence_pipeline.md
  Model1_Methods_copy_ready.md
  assets/model1_pipeline_with_wgs.svg

docs/
  WGS_filtering.md
  installation, troubleshooting, outputs and limitations

results/
  placeholders only; no numerical result files are committed yet
```

## RNA-seq quick start

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
STAR_INDEX=/data/ydx/igem/STAR_index
REDITOOLS=/data/ydx/igem/REDItools2

python3 scripts/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest config/samples.tsv \
  --threads 16

python3 scripts/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest config/samples.tsv \
  --threads 50

python3 scripts/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest config/samples.tsv \
  --java-options=-Xmx16g

conda activate reditools2_py2
nohup bash scripts/run_reditools_all_samples.sh \
  "$PROJECT" "$REF" "$REDITOOLS" "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

## Public HEK293T WGS quick start

Configured runs:

```text
SRR37832939
SRR37832940
SRR37832941
```

Create the environment and inspect metadata:

```bash
conda env create -f environment/wgs_pipeline.yml
conda activate rewire_wgs

bash scripts/wgs/00_check_sra_metadata.sh \
  config/wgs_runs.tsv \
  wgs_sra_metadata.tsv
```

Run automatically:

```bash
WGS_OUT=/data/ydx/igem/HEK293T_public_WGS_3runs

nohup bash scripts/wgs/run_3run_wgs_pipeline.sh \
  --runs config/wgs_runs.tsv \
  --reference "$REF" \
  --outdir "$WGS_OUT" \
  --mode auto \
  --threads 32 \
  --min-dp 10 \
  --min-alt 3 \
  --min-vaf 0.05 \
  --min-qual 20 \
  > "$WGS_OUT.pipeline.log" 2>&1 &
```

The metadata determine whether runs are merged as one BioSample or called separately to produce a two-of-three exact-allele consensus blacklist.

## Final integration

```bash
python3 scripts/reditools_union_to_vcf.py \
  --manifest config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"

python3 scripts/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache /path/to/vep_cache

bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  config/samples.tsv

python3 scripts/filter_c_to_u_and_compare.py \
  --manifest config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --wgs-vcf "$WGS_OUT/vcf/HEK293T_3runs.consensus2of3.SNV.vcf.gz" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

If metadata resolve to merge mode, use:

```text
$WGS_OUT/vcf/HEK293T_3runs.filtered.SNV.vcf.gz
```

## Key methodological decisions

### REDItools2 discovery

```text
-S       report positions containing observed substitutions
-me 20   require at least 20 edited reads at a reported site
-q 20    minimum mapping quality
-bq 30   minimum base quality
```

### Transcript-oriented interpretation

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

### Conservative evidence rule

```text
3/3 treated calls
0/3 control calls
candidate-site depth ≥20 in all six samples
strand consistent with C-to-U
exact allele absent from the selected public WGS blacklist
```

### WGS single-run thresholds

```text
DP ≥10
ALT reads ≥3
VAF ≥0.05
QUAL ≥20
```

### Public WGS boundary

The WGS data are external public HEK293T genomes. They strengthen genomic-SNV filtering but cannot be described as matched WGS from the exact experimental cell batch.

## Results status

Methods, scripts, environments and editable Wiki figures are included. Numerical result files and result-specific plots are intentionally excluded until the six RNA-seq samples and selected WGS workflow pass the same integrity checks.

The repository state before the Model 1 rebuild remains available on branch `backup-before-paper-pipeline-20260704`.
