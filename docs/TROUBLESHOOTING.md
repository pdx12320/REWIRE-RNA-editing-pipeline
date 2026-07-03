# Troubleshooting

## Missing MPI slots

Open MPI may expose fewer physical-core slots than the number of logical CPUs shown by `lscpu`. Use the allocation reported by Open MPI. On the current server, 32 slots are available, so 30 REDItools2 processes are used.

## Missing Python 2 modules

REDItools2 must be launched with the Python interpreter from its dedicated environment. That interpreter must import mpi4py, pysam, sortedcontainers, psutil, and netifaces. The mpi4py package must be compiled against the same MPI installation used at runtime.

## Coverage file is null

The parallel REDItools2 program requires both a complete coverage file and the directory containing per-contig coverage files. The pipeline passes these through the `-G` and `-D` arguments.

## Long GL and KI coverage logs

This is expected. GRCh38 includes supplementary contigs with GL and KI identifiers. The coverage step iterates through every entry in the FASTA index.

## ncurses warnings from samtools

Warnings about missing symbol-version information in `libncursesw` or `libtinfow` usually mean that Conda libraries are shadowing system libraries. They are non-fatal when samtools continues writing output. A known system samtools binary can be used to avoid the warning.

## Missing read group

A GATK null-pointer error involving `getReadGroup()` indicates a BAM without read-group tags. Add read groups before MarkDuplicates or rerun STAR with read-group metadata.

## Basic validation

Use `samtools quickcheck` on every final SplitNCigarReads BAM. Confirm that each REDItools2 table is non-empty, can be decompressed, starts with the expected header, and has a matching tabix index.
