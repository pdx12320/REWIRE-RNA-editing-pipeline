# Pipeline implementation

This directory contains the executable implementation behind the copy-ready iGEM text in [`../wiki/README.md`](../wiki/README.md).

## Contents

```text
pipeline/
├── README.md
├── OUTPUTS.md
├── TROUBLESHOOTING.md
├── config/
│   ├── samples.tsv
│   └── wgs_runs.tsv
├── env/
│   ├── reditools2_py2.yml
│   └── wgs_pipeline.yml
└── scripts/
    ├── rna/
    │   ├── download_sra_fastq.py
    │   ├── run_star_alignment.py
    │   ├── run_gatk_preprocessing.py
    │   ├── generate_reditools_coverage_limited.sh
    │   ├── run_reditools_all_samples.sh
    │   ├── rebuild_reditools_file_list.py
    │   ├── reditools_union_to_vcf.py
    │   ├── run_vep_annotation.py
    │   ├── build_candidate_depth_tables.sh
    │   ├── filter_c_to_u_and_compare.py
    │   ├── filter_calls.py
    │   └── filter_utils.py
    └── wgs/
        ├── 00_check_sra_metadata.sh
        ├── run_3run_wgs_pipeline.sh
        ├── filter_single_sample_vcf.py
        └── build_consensus_blacklist.py
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

### WGS branch

```bash
conda env create -f pipeline/env/wgs_pipeline.yml
conda activate rewire_wgs
```

### REDItools2 environment

The provided YAML creates a minimal Python 2.7 environment:

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

## 2. Fixed inputs

RNA samples are stored in [`config/samples.tsv`](config/samples.tsv). Public HEK293T WGS runs are stored in [`config/wgs_runs.tsv`](config/wgs_runs.tsv).

All branches must use the same GRCh38 FASTA and chromosome naming convention. Do not mix `chr1` and `1`, or GRCh37 and GRCh38 coordinates, without explicit normalisation.

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

## 4. Public HEK293T WGS workflow

Check whether the three runs share one BioSample:

```bash
bash pipeline/scripts/wgs/00_check_sra_metadata.sh \
  pipeline/config/wgs_runs.tsv \
  wgs_sra_metadata.tsv
```

Run the complete workflow:

```bash
conda activate rewire_wgs
WGS_OUT=/data/ydx/igem/HEK293T_public_WGS_3runs

nohup bash pipeline/scripts/wgs/run_3run_wgs_pipeline.sh \
  --runs pipeline/config/wgs_runs.tsv \
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

`auto` merges runs only when the metadata indicate one BioSample. Otherwise it creates per-run call sets and a two-of-three exact-allele consensus blacklist.

## 5. Final integration

Choose the WGS blacklist produced by the resolved mode:

```text
merge mode:
  $WGS_OUT/vcf/HEK293T_3runs.filtered.SNV.vcf.gz

consensus mode:
  $WGS_OUT/vcf/HEK293T_3runs.consensus2of3.SNV.vcf.gz
```

Run the final filter:

```bash
python3 pipeline/scripts/rna/filter_c_to_u_and_compare.py \
  --manifest pipeline/config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --wgs-vcf "$WGS_BLACKLIST" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

## 6. Reproducibility rules

- Record software versions and reference checksums.
- Use one GRCh38 FASTA across STAR, GATK, REDItools2, BWA and bcftools.
- Preserve complete supplementary-contig identifiers, including `.1` and `.2` suffixes.
- Do not interpret an absent control call without candidate-site depth.
- Match WGS evidence by exact `CHROM:POS:REF:ALT`, not coordinate alone.
- Keep raw data and large intermediate files outside GitHub.

See [`OUTPUTS.md`](OUTPUTS.md) for expected files and [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md) for known deployment issues.
