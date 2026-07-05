# Pipeline implementation

This directory contains the executable implementation behind the iGEM-ready explanation in [`../wiki/README.md`](../wiki/README.md).

## One-sentence workflow

Six RNA-seq libraries generate replicate-, control-, depth- and strand-aware C-to-U evidence; a database-released 293T variant catalogue is converted from hg18 to GRCh38 and joined to the RNA evidence by exact allele.

## Contents

```text
pipeline/
├── README.md
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
    │   └── process_293T_CG_to_GRCh38.sh
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

## 1. Required software

### RNA-seq branch

- SRA Toolkit
- STAR
- samtools and htslib
- GATK 4 with Java 17
- REDItools2
- Python 2.7 with `mpi4py`, `pysam`, `sortedcontainers`, `psutil` and `netifaces`
- Open MPI
- VEP with an offline GRCh38 cache
- Python 3

### 293T catalogue branch

```bash
conda env create -f pipeline/env/genomic_catalogue.yml
conda activate rewire_catalogue
```

This environment provides CrossMap, bcftools, samtools, bgzip and tabix. BWA, GATK and SRA Toolkit are not required for the catalogue branch because the input is an existing VCF rather than raw WGS reads.

### REDItools2 environment

```bash
conda env create -f pipeline/env/reditools2_py2.yml
conda activate reditools2_py2
```

If the Python-2 package solver cannot install the remaining modules, bootstrap the Python-2 pip installer and use versions compatible with Python 2.7:

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

Use the same MPI implementation for `mpirun`, `mpicc` and the compiled `mpi4py` module.

## 2. Fixed inputs and reference rule

RNA samples are stored in [`config/samples.tsv`](config/samples.tsv).

The genomic catalogue input is `293T_CG.vcf.gz` from the [HEK293 Genome Project data page](https://hek293genome.org/v2/data.php). The source file uses NCBI build 36/hg18 and must be converted before comparison with GRCh38 RNA-seq results.

All modules must use the same GRCh38 FASTA and chromosome naming convention. Do not mix `chr1` and `1`, or hg18 and GRCh38 coordinates, without explicit harmonization.

## 3. RNA-seq workflow

Set paths:

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
STAR_INDEX=/data/ydx/igem/STAR_index
REDITOOLS=/data/ydx/igem/REDItools2
```

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

### GATK RNA preprocessing

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

The three trailing values are MPI processes, concurrent coverage jobs and compression threads.

### Build union candidate coordinates

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

### Candidate depth in all samples

```bash
bash pipeline/scripts/rna/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  pipeline/config/samples.tsv
```

## 4. HEK293 Genome Project catalogue workflow

Download the following two files outside the repository:

```text
293T_CG.vcf.gz                  HEK293 Genome Project data page
hg18ToHg38.over.chain.gz        UCSC liftOver chain
```

The processing script accepts either a text VCF or a gzip-compressed VCF, including a compressed file that was saved with a `.vcf` extension.

```bash
conda activate rewire_catalogue

CATALOGUE_IN=/data/ydx/igem/293T_CG.vcf
CHAIN=/data/ydx/igem/hg18ToHg38.over.chain.gz
CATALOGUE_OUT=/data/ydx/igem/293T_CG_GRCh38

nohup bash pipeline/scripts/catalogue/process_293T_CG_to_GRCh38.sh \
  --input "$CATALOGUE_IN" \
  --reference "$REF" \
  --chain "$CHAIN" \
  --outdir "$CATALOGUE_OUT" \
  --threads 16 \
  > "$CATALOGUE_OUT.nohup.log" 2>&1 &
```

The script performs:

```text
input-format normalization
→ PASS biallelic SNP selection
→ CrossMap hg18-to-GRCh38 liftover
→ GRCh38 REF validation and mismatch removal
→ normalization and coordinate sorting
→ bgzip/tabix indexing
→ C>T / G>A catalogue extraction
→ QC summary and checksums
```

Recommended integration file:

```text
$CATALOGUE_OUT/293T_CG.GRCh38.CtoU_relevant.SNV.vcf.gz
```

See [`CATALOGUE_PROVENANCE.md`](CATALOGUE_PROVENANCE.md) for the observed QC counts and interpretation boundary.

## 5. Final integration

Run the final filter:

```bash
CATALOGUE_VCF="$CATALOGUE_OUT/293T_CG.GRCh38.CtoU_relevant.SNV.vcf.gz"

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

The older `--wgs-vcf` argument remains as a compatibility alias, but new analyses should use `--variant-catalogue-vcf`.

## 6. Reproducibility rules

- Record software versions, source URLs and file checksums.
- Use one GRCh38 FASTA across STAR, GATK, REDItools2, CrossMap output validation and bcftools normalization.
- Preserve complete supplementary-contig identifiers, including `.1` and `.2` suffixes.
- Do not interpret an absent control call without candidate-site depth.
- Match catalogue evidence by exact `CHROM:POS:REF:ALT`, not coordinate alone.
- Keep catalogue-overlapping records in the complete site matrix even when they are excluded from the high-confidence set.
- Keep raw data and large intermediate files outside GitHub.

See [`OUTPUTS.md`](OUTPUTS.md) for expected files and [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) for known deployment issues.
