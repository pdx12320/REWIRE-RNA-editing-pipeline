# Troubleshooting

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

Do not treat the site as control-negative until candidate-site depth has been queried from the BAM. REDItools2 `-me 20` can suppress a site even when low-level alternate reads are present.

## `bcftools` reports `Exec format error` for `293T_CG.vcf`

The file may be gzip-compressed even though it has a `.vcf` extension, or may contain bytes before the gzip header. The catalogue script detects the content format and writes a normalized temporary VCF before processing. Do not infer compression from the filename alone.

Check manually with:

```bash
file 293T_CG.vcf
head -c 16 293T_CG.vcf | od -An -t x1
```

A gzip stream begins with `1f 8b 08`; a text VCF begins with `##fileformat`.

## UCSC chain download times out

Download `hg18ToHg38.over.chain.gz` on another machine and copy it to the server. Verify it before rerunning:

```bash
gzip -t /path/hg18ToHg38.over.chain.gz
```

Pass the local file explicitly with `--chain`.

## CrossMap output cannot be indexed because positions are unsorted

CrossMap does not guarantee final VCF coordinate order. The catalogue script runs `bcftools sort` after REF validation and before tabix indexing. Do not index the unsorted CrossMap output directly.

Typical error:

```text
Unsorted positions on sequence ...
index: failed to create index
```

## Liftover produces REF mismatches

Some mapped records may not match the target GRCh38 REF allele. The script uses:

```bash
bcftools norm -f GRCh38.fa -c x
```

to remove incompatible records. The frozen project conversion removed 22,761 REF-mismatch records. A substantially different count should trigger a reference and chain audit.

## Catalogue variant count is unexpectedly small

The final normalized catalogue should contain millions of SNPs, not hundreds. Check:

```bash
bcftools index -n 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
```

The frozen conversion retained 2,885,725 SNPs. A small count usually indicates the wrong input file, an assembly mismatch or an inappropriate genotype/filter expression.
