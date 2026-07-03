# Model 1 — Code and Reproducibility

This page collects the executable commands behind the RNA-editing evidence pipeline. It is written for an iGEM Software/Dry Lab page: each section explains what the command does, what it produces, and what must be checked before moving on.

---

## 1. Define the project paths

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
STAR_INDEX=/data/ydx/igem/STAR_index
REDITOOLS=/data/ydx/igem/REDItools2
MANIFEST=config/samples.tsv
```

Reference-dependent files must use the same GRCh38 naming convention:

```text
reference FASTA
reference FASTA index (.fai)
GATK sequence dictionary (.dict)
STAR index
VEP cache
optional WGS VCF
```

---

## 2. Download SRA data and create FASTQ files

```bash
python3 scripts/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest "$MANIFEST" \
  --threads 16
```

Expected files:

```text
fastq/CU517_GC_T1_1.fastq.gz
fastq/CU517_GC_T1_2.fastq.gz
...
fastq/CU517_GC_C3_2.fastq.gz
```

Checkpoint:

```bash
find "$PROJECT/fastq" -name '*.fastq.gz' -size +0c | wc -l
```

Expected count: 12 non-empty FASTQ files.

---

## 3. STAR alignment

```bash
python3 scripts/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest "$MANIFEST" \
  --threads 50
```

Checkpoint:

```bash
for bam in "$PROJECT"/bam/star/*.bam; do
  samtools quickcheck -v "$bam" || echo "FAILED: $bam"
done
```

Read-group check:

```bash
samtools view -H "$PROJECT/bam/star/CU517_GC_T1.Aligned.sortedByCoord.out.bam" \
  | grep '^@RG'
```

---

## 4. GATK RNA preprocessing

```bash
python3 scripts/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest "$MANIFEST" \
  --java-options=-Xmx16g
```

Expected final BAM files:

```text
bam/splitncigarreads/CU517_GC_T1.splitncigarreads.bam
...
bam/splitncigarreads/CU517_GC_C3.splitncigarreads.bam
```

Checkpoint:

```bash
for bam in "$PROJECT"/bam/splitncigarreads/*.bam; do
  samtools quickcheck -v "$bam" || echo "FAILED: $bam"
done
```

---

## 5. REDItools2 Python environment

REDItools2 uses a Python 2 environment. Required modules include:

```text
mpi4py
pysam
sortedcontainers
psutil
netifaces
```

Check the interpreter:

```bash
conda activate reditools2_py2

"$CONDA_PREFIX/bin/python" - <<'PY'
from mpi4py import MPI
import pysam, sortedcontainers, psutil, netifaces
print("REDItools2 environment OK")
print(MPI.Get_library_version())
PY
```

Check MPI:

```bash
mpirun -np 2 "$CONDA_PREFIX/bin/python" -c \
'from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank(), MPI.COMM_WORLD.Get_size())'
```

---

## 6. Generate coverage with limited concurrency

The helper script creates one coverage file per contig and one complete coverage file per sample.

```bash
bash scripts/generate_reditools_coverage_limited.sh \
  "$PROJECT/bam/splitncigarreads/CU517_GC_T1.splitncigarreads.bam" \
  "$PROJECT/reditools/coverage/CU517_GC_T1" \
  "${REF}.fai" \
  8
```

Monitor progress:

```bash
pgrep -af 'samtools depth'
du -sh "$PROJECT/reditools/coverage/CU517_GC_T1"
```

The appearance of GL and KI contigs is expected. Their `.1` or `.2` suffixes must remain unchanged.

---

## 7. Run REDItools2 for all six samples

```bash
nohup bash scripts/run_reditools_all_samples.sh \
  "$PROJECT" \
  "$REF" \
  "$REDITOOLS" \
  "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools_all_samples.log" 2>&1 &
```

The final three numbers mean:

```text
30  MPI processes
8   concurrent coverage jobs
8   compression threads
```

Monitor the global log:

```bash
tail -f "$PROJECT/logs/reditools_all_samples.log"
```

Monitor one sample:

```bash
grep -E 'COVERAGE RECEIVED|RECEIVED IM_FREE|WHOLE PARALLEL|ERROR' \
  "$PROJECT/logs/CU517_GC_T1.REDItools2.log" | tail -30
```

Expected per-sample files:

```text
reditools/tables/CU517_GC_T1.txt.gz
reditools/tables/CU517_GC_T1.txt.gz.tbi
```

Integrity checks:

```bash
gzip -t "$PROJECT/reditools/tables/CU517_GC_T1.txt.gz"

tabix -l "$PROJECT/reditools/tables/CU517_GC_T1.txt.gz" | head

zcat "$PROJECT/reditools/tables/CU517_GC_T1.txt.gz" | head -n 3
```

---

## 8. REDItools2 contig-name patch

The original temporary-file parser removed the contig version suffix. We replaced:

```python
pieces = re.sub("\\..*", "", os.path.basename(little_file)).split("#")
```

with:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

Check the installed source:

```bash
grep -n 'pieces = ' \
  /data/ydx/igem/REDItools2/src/cineca/parallel_reditools.py | tail
```

Expected code:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

---

## 9. Build the union VCF and candidate BED

```bash
mkdir -p "$PROJECT/vcf"

python3 scripts/reditools_union_to_vcf.py \
  --manifest "$MANIFEST" \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"
```

The VCF is used for strand annotation. The BED file is used to query candidate-site depth in every sample.

---

## 10. VEP transcript-strand annotation

```bash
python3 scripts/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache /path/to/vep_cache
```

Expected output:

```text
vep/CU5.17_EGFP_GC.vep.tsv
```

Required fields:

```text
Uploaded_variation
Location
Allele
STRAND
```

---

## 11. Measure candidate-site depth in all replicates

```bash
bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  "$MANIFEST"
```

The script uses:

```bash
samtools depth -a -q 30 -Q 20 -b union_candidates.bed sample.bam
```

In `samtools depth`:

```text
-q 30  minimum base quality
-Q 20  minimum mapping quality
-a     include requested positions with zero depth
```

---

## 12. Build the final evidence matrix

```bash
python3 scripts/filter_c_to_u_and_compare.py \
  --manifest "$MANIFEST" \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

Optional matched WGS filtering:

```bash
--wgs-vcf /path/to/HEK293T.filtered.vcf.gz
```

No result tables are committed to this wiki package until all six samples complete the workflow.

---

## Repository code map

| Stage | Script |
|---|---|
| SRA download | `scripts/download_sra_fastq.py` |
| Alignment | `scripts/run_star_alignment.py` |
| GATK preprocessing | `scripts/run_gatk_preprocessing.py` |
| Coverage generation | `scripts/generate_reditools_coverage_limited.sh` |
| REDItools2 | `scripts/run_reditools_all_samples.sh` |
| Union variants | `scripts/reditools_union_to_vcf.py` |
| VEP | `scripts/run_vep_annotation.py` |
| All-sample depth | `scripts/build_candidate_depth_tables.sh` |
| C-to-U filtering | `scripts/filter_utils.py`, `scripts/filter_calls.py` |
| Evidence matrix | `scripts/filter_c_to_u_and_compare.py` |
