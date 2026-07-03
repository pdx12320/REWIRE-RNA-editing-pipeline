# Model 1 code and commands

This page contains the concrete commands used by the RNA-editing evidence pipeline. It can be linked from the iGEM Wiki as the reproducibility section.

## 1. Define paths

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
STAR_INDEX=/data/ydx/igem/STAR_index
REDITOOLS=/data/ydx/igem/REDItools2
VEP_CACHE=/path/to/vep_cache
MANIFEST=config/samples.tsv
```

## 2. Download SRA and convert to FASTQ

```bash
python3 scripts/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest "$MANIFEST" \
  --threads 16
```

## 3. STAR two-pass alignment

```bash
python3 scripts/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest "$MANIFEST" \
  --threads 50
```

Expected BAM pattern:

```text
$PROJECT/bam/star/CU517_GC_T1.Aligned.sortedByCoord.out.bam
```

## 4. GATK RNA preprocessing

```bash
python3 scripts/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest "$MANIFEST" \
  --java-options=-Xmx16g
```

Expected final preprocessing BAM pattern:

```text
$PROJECT/bam/splitncigarreads/CU517_GC_T1.splitncigarreads.bam
```

## 5. Activate the REDItools2 environment

```bash
conda activate reditools2_py2

python - <<'PY'
from mpi4py import MPI
import pysam, sortedcontainers, psutil, netifaces
print("REDItools2 dependencies: OK")
print(MPI.Get_library_version())
PY
```

## 6. Run REDItools2 for all six samples

```bash
nohup bash scripts/run_reditools_all_samples.sh \
  "$PROJECT" \
  "$REF" \
  "$REDITOOLS" \
  "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

Parameter interpretation:

```text
30 = MPI processes
8  = simultaneous samtools depth jobs
8  = bgzip compression threads
```

Monitor the wrapper log:

```bash
tail -f "$PROJECT/logs/reditools.log"
```

Monitor one sample:

```bash
grep -E 'COVERAGE RECEIVED|RECEIVED IM_FREE|WHOLE PARALLEL|ERROR' \
  "$PROJECT/logs/CU517_GC_T1.REDItools2.log" | tail -30
```

Check active processes:

```bash
pgrep -af 'samtools depth|parallel_reditools.py|mpirun'
```

## 7. Verify each REDItools2 table

```bash
for table in "$PROJECT"/reditools/tables/CU517_GC_*.txt.gz; do
  echo "Checking $table"
  gzip -t "$table"
  zcat "$table" | awk -F '\t' '
    NR==1 && NF!=14 {print "bad header", NF; exit 1}
    NR>1 && NF!=14 {print "bad row", NR, NF; exit 1}
    END {print "rows", NR-1}
  '
done
```

Check tabix indexes:

```bash
for table in "$PROJECT"/reditools/tables/CU517_GC_*.txt.gz; do
  test -s "${table}.tbi" || echo "Missing index: ${table}.tbi"
done
```

## 8. Build the union VCF and BED

```bash
mkdir -p "$PROJECT/vcf"

python3 scripts/reditools_union_to_vcf.py \
  --manifest "$MANIFEST" \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"
```

## 9. Run offline VEP strand annotation

```bash
python3 scripts/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache "$VEP_CACHE"
```

Expected output:

```text
$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv
```

## 10. Measure candidate depth in all six samples

```bash
bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  "$MANIFEST"
```

Depth is measured with:

```text
base quality ≥30
mapping quality ≥20
```

## 11. Build the final evidence matrix

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

With matched HEK293T WGS filtering:

```bash
python3 scripts/filter_c_to_u_and_compare.py \
  --manifest "$MANIFEST" \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --wgs-vcf /path/to/HEK293T.filtered.vcf.gz \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

## 12. Expected output structure

```text
$PROJECT/
├── fastq/
├── bam/
│   ├── star/
│   ├── readgroups/
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

The repository does not currently include numerical result files. Only source code, documentation, and Wiki-ready text are committed.
