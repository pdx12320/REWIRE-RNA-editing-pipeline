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

Warnings about `libncursesw.so.6` or `libtinfow.so.6` usually reflect Conda library shadowing. If `samtools` continues to produce valid output, the warning is not itself a pipeline failure. A system samtools build can avoid the warning.

## A control has no REDItools2 call

Do not treat the site as control-negative until candidate-site depth has been queried from the BAM. REDItools2 `-me 20` can suppress a site even when low-level alternate reads are present.

## WGS runs are consecutive accessions

Consecutive SRR identifiers do not prove that runs belong to one biological sample. Run `00_check_sra_metadata.sh`; merge only when the BioSample accessions agree.
