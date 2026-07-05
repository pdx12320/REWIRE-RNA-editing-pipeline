# Troubleshooting

These fixes correspond to observed failures during development. The related DBTL decisions are recorded in [`../docs/ENGINEERING_CYCLE.md`](../docs/ENGINEERING_CYCLE.md).

## GATK reports a null read group

Symptom:

```text
SAMRecord.getReadGroup() is null
```

Confirm that the BAM header contains `@RG`. Add or replace read groups before `MarkDuplicates` when they are absent.

## REDItools2 cannot import Python modules

Confirm that the Python interpreter passed to REDItools2 contains:

```text
mpi4py
pysam
sortedcontainers
psutil
netifaces
```

Pass the interpreter explicitly:

```bash
"$CONDA_PREFIX/bin/python"
```

## Open MPI reports insufficient slots

Inspect the allocation:

```bash
mpirun --display-allocation -np 1 hostname
```

Do not request more MPI processes than the available slots. On a 32-slot allocation, 30 workers leave headroom for the controller and system processes.

## Coverage generation appears to stop at GL/KI contigs

GL and KI identifiers are supplementary GRCh38 contigs. Their presence is expected. Check active `samtools depth` processes and directory growth before assuming the job has stalled.

## REDItools2 fails on `chrGL000009`

Cause: version suffixes such as `.1` and `.2` were removed from temporary filenames.

Required parser:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

This removes only `.gz` and preserves the complete contig identifier.

## samtools reports ncurses version warnings

Warnings about `libncursesw.so.6` or `libtinfow.so.6` usually reflect Conda library shadowing. If `samtools` continues to produce valid output, the warning is not itself a pipeline failure. A dedicated environment or system samtools build avoids the warning.

## A control has no REDItools2 call

Do not treat the site as control-negative until candidate-site depth and base counts have been queried from the BAM. REDItools2 `-me 20` can suppress a site even when lower-level alternate reads are present.

For the frozen legacy result, control non-calls are retained as a stated limitation. The 3,333 retained rows are screening candidates, not a fully depth-qualified set.

## `candidate_depth/` is missing

The strict integration route requires:

```text
$PROJECT/candidate_depth/<sample>.candidate_depth.tsv.gz
```

Generate these files with:

```bash
bash pipeline/scripts/rna/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  pipeline/config/samples.tsv
```

If a completed legacy `treatment_specific.tsv.gz` already exists, use the compatibility route instead:

```bash
python3 pipeline/scripts/catalogue/filter_existing_treatment_specific_by_293T.py \
  --treatment-specific "$PROJECT/final/CU5.17_EGFP_GC.treatment_specific.tsv.gz" \
  --catalogue-vcf /path/293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz \
  --output-dir "$PROJECT/final_with_293T_catalogue"
```

This route does not fabricate missing depth evidence.

## The site matrix lacks `all_replicates_depth_pass`

The table was produced by an older helper script. Do not rename another column or fill the missing value manually.

Use `filter_existing_treatment_specific_by_293T.py`, or regenerate candidate-depth tables and rerun the strict integration route from the original REDItools and VEP inputs.

## `bcftools` reports `Exec format error` for `293T_CG.vcf`

The file may be gzip-compressed even though it has a `.vcf` extension. The catalogue script detects the content format before processing.

Check manually with:

```bash
file 293T_CG.vcf
head -c 16 293T_CG.vcf | od -An -t x1
```

A gzip stream begins with `1f 8b 08`; a text VCF begins with `##fileformat`.

## UCSC chain download times out

Download `hg18ToHg38.over.chain.gz` on another machine, copy it to the server and verify it:

```bash
gzip -t /path/hg18ToHg38.over.chain.gz
```

Pass the local file explicitly with `--chain`.

## CrossMap output cannot be indexed because positions are unsorted

CrossMap does not guarantee coordinate order. The catalogue script runs `bcftools sort` before tabix indexing.

Typical error:

```text
Unsorted positions on sequence ...
index: failed to create index
```

## Liftover produces REF mismatches

Some mapped records do not match the target GRCh38 REF allele. The script removes incompatible records with:

```bash
bcftools norm -f GRCh38.fa -c x
```

The frozen conversion removed 22,761 REF-mismatch records.

## Catalogue variant count is unexpectedly small

The final catalogue should contain millions of SNPs:

```bash
bcftools index -n 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
```

The frozen conversion retained 2,885,725 SNPs. A much smaller count usually indicates the wrong input file, an assembly mismatch or an inappropriate filter.
