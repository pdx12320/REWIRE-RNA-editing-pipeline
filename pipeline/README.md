# Pipeline implementation

This document gives the executable workflow. The reasoning, failed tests and DBTL decisions are recorded in [`../docs/ENGINEERING_CYCLE.md`](../docs/ENGINEERING_CYCLE.md).

## Data flow

```text
six RNA-seq libraries
→ GRCh38 + sequence-verified EGFP-GC reporter reference
→ uniform STAR and GATK preprocessing for all six libraries
→ REDItools2 substitution calls
→ VEP transcript orientation
→ treated/control comparison
→ exact-allele filtering against the GRCh38 293T catalogue
```

## Strict CU5.17 model-training rebuild

For the corrected training dataset, run the end-to-end strict workflow:

```bash
bash pipeline/scripts/rna/rebuild_strict_cu517_dataset.sh \
  "$PROJECT" "$GRCH38_FASTA" "$GENCODE_GTF" \
  "$SEQUENCE_VERIFIED_EGFP_GC_REPORTER_FASTA" \
  pipeline/config/samples.tsv 32
```

This workflow verifies EGFP positions 458–459 as `GC`, builds an augmented STAR
reference, recreates all six BAMs through one preprocessing route, requires
`NH=1`, audits sequence mappability, enumerates every covered exonic cytidine
with a contiguous exon-contained 101-nt window,
keeps depth `>50` in every replicate, constructs zero-alt strict negatives and
positive sites, then writes a 1:200 gene-disjoint handoff. See
[`LAMAR_TRAINING_LABELS.md`](LAMAR_TRAINING_LABELS.md) for definitions and output
files.

The exact reporter FASTA must come from the sequenced plasmid/construct or the
authors. The paper and its one-file public analysis repository do not provide a
complete reporter sequence, so this repository does not substitute a guessed
pEGFP-C1 sequence.

## Repository structure

```text
pipeline/
├── CATALOGUE_PROVENANCE.md
├── OUTPUTS.md
├── TROUBLESHOOTING.md
├── config/
│   └── samples.tsv
├── env/
│   ├── reditools2_py2.yml
│   └── genomic_catalogue.yml
└── scripts/
    ├── catalogue/
    │   ├── process_293T_CG_to_GRCh38.sh
    │   └── filter_existing_treatment_specific_by_293T.py
    └── rna/
        ├── download_sra_fastq.py
        ├── run_star_alignment.py
        ├── run_gatk_preprocessing.py
        ├── generate_reditools_coverage_limited.sh
        ├── run_reditools_all_samples.sh
        ├── rebuild_reditools_file_list.py
        ├── reditools_union_to_vcf.py
        ├── run_vep_annotation.py
        ├── build_candidate_depth_tables.sh
        ├── filter_c_to_u_and_compare.py
        ├── filter_calls.py
        └── filter_utils.py
```

## 1. Software

### RNA branch

- SRA Toolkit
- STAR
- samtools and htslib
- BWA 0.7.17 or newer (101-nt mappability audit)
- GATK 4 with Java 17
- REDItools2
- Open MPI
- Python 2.7 with `mpi4py`, `pysam`, `sortedcontainers`, `psutil` and `netifaces`
- VEP with an offline GRCh38 cache
- Python 3

### Catalogue branch

```bash
conda env create -f pipeline/env/genomic_catalogue.yml
conda activate rewire_catalogue
```

This environment provides CrossMap, bcftools, samtools, bgzip and tabix.

### REDItools2 environment

```bash
conda env create -f pipeline/env/reditools2_py2.yml
conda activate reditools2_py2
```

When the Python-2 solver cannot install the remaining modules:

```bash
curl -L https://bootstrap.pypa.io/pip/2.7/get-pip.py -o /tmp/get-pip-py2.py
python /tmp/get-pip-py2.py
python -m pip install "setuptools<45" "wheel<0.35"
python -m pip install \
  "pysam==0.15.4" \
  "sortedcontainers==2.2.2" \
  "psutil==5.6.7" \
  "netifaces==0.10.9"

export MPICC="$(command -v mpicc)"
python -m pip install --no-cache-dir --no-binary=mpi4py "mpi4py==3.0.3"
```

Use the same MPI implementation for `mpirun`, `mpicc` and `mpi4py`.

## 2. Fixed inputs

RNA samples are listed in [`config/samples.tsv`](config/samples.tsv).

Set paths:

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
STAR_INDEX=/data/ydx/igem/STAR_index
REDITOOLS=/data/ydx/igem/REDItools2
```

Use the same GRCh38 FASTA and chromosome naming convention across all modules.

## 3. RNA-seq workflow

### Download FASTQ

```bash
python3 pipeline/scripts/rna/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest pipeline/config/samples.tsv \
  --threads 16
```

### STAR alignment

```bash
python3 pipeline/scripts/rna/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest pipeline/config/samples.tsv \
  --threads 50
```

### GATK preprocessing

```bash
python3 pipeline/scripts/rna/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest pipeline/config/samples.tsv \
  --java-options=-Xmx16g
```

### REDItools2

```bash
conda activate reditools2_py2

nohup bash pipeline/scripts/rna/run_reditools_all_samples.sh \
  "$PROJECT" "$REF" "$REDITOOLS" "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

The final three values are MPI processes, concurrent coverage jobs and compression threads.

### Union candidate coordinates

```bash
python3 pipeline/scripts/rna/reditools_union_to_vcf.py \
  --manifest pipeline/config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"
```

### VEP annotation

```bash
python3 pipeline/scripts/rna/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache /path/to/vep_cache
```

### Strict candidate-depth measurement

```bash
bash pipeline/scripts/rna/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  pipeline/config/samples.tsv
```

This step measures each candidate directly in all six BAM files. It is required for the strict integration route.

## 4. 293T catalogue workflow

Input files:

```text
293T_CG.vcf.gz                  HEK293 Genome Project
hg18ToHg38.over.chain.gz        UCSC liftOver chain
```

The source VCF uses NCBI build 36/hg18 and must be converted before comparison with GRCh38 RNA coordinates.

```bash
conda activate rewire_catalogue

CATALOGUE_IN=/data/ydx/igem/293T_CG.vcf
CHAIN=/data/ydx/igem/hg18ToHg38.over.chain.gz
CATALOGUE_OUT=/data/ydx/igem/293T_CG_GRCh38

bash pipeline/scripts/catalogue/process_293T_CG_to_GRCh38.sh \
  --input "$CATALOGUE_IN" \
  --reference "$REF" \
  --chain "$CHAIN" \
  --outdir "$CATALOGUE_OUT" \
  --threads 16
```

The conversion performs:

```text
content-based input detection
→ PASS biallelic SNP selection
→ CrossMap hg18-to-GRCh38 liftover
→ GRCh38 REF validation
→ normalization and coordinate sorting
→ bgzip and tabix indexing
→ QC summary and checksums
```

Frozen QC:

```text
source PASS biallelic SNPs       2,914,465
CrossMap-unmapped records            5,979
GRCh38 REF mismatches removed       22,761
final GRCh38 SNPs                2,885,725
```

See [`CATALOGUE_PROVENANCE.md`](CATALOGUE_PROVENANCE.md).

## 5. Final integration

### Route A — strict integration

Use this route when all six candidate-depth tables exist:

```bash
CATALOGUE_VCF="$CATALOGUE_OUT/293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz"

python3 pipeline/scripts/rna/filter_c_to_u_and_compare.py \
  --manifest pipeline/config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --variant-catalogue-vcf "$CATALOGUE_VCF" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

### Route B — completed legacy table

Use this route when a completed treatment-specific table exists but the older site matrix lacks `all_replicates_depth_pass` or the `candidate_depth/` directory is unavailable:

```bash
INPUT="$PROJECT/final/CU5.17_EGFP_GC.treatment_specific.tsv.gz"
CATALOGUE_VCF="$CATALOGUE_OUT/293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz"
OUT="$PROJECT/final_with_293T_catalogue"

python3 pipeline/scripts/catalogue/filter_existing_treatment_specific_by_293T.py \
  --treatment-specific "$INPUT" \
  --catalogue-vcf "$CATALOGUE_VCF" \
  --output-dir "$OUT"
```

The compatibility route adds `genomic_catalogue_overlap`, preserves excluded alleles and does not invent missing depth evidence.

Frozen result:

```text
treatment-specific before catalogue    3,349
exact catalogue overlaps                  16
retained screening candidates          3,333
```

## 6. Reproducibility rules

1. Use one GRCh38 FASTA across the complete workflow.
2. Record software versions, input URLs and checksums.
3. Preserve supplementary-contig version suffixes.
4. Do not equate a missing control call with zero editing.
5. Match the catalogue by exact `CHROM:POS:REF:ALT`.
6. Keep catalogue-overlapping records in a separate exclusion table.
7. Do not describe the legacy retained set as fully depth-qualified.

See [`OUTPUTS.md`](OUTPUTS.md) and [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md).
