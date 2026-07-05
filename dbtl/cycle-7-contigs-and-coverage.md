# Cycle 7 — Preserve supplementary contigs and coverage output

## Design

### Question

Could coverage generation and temporary-file parsing handle the complete GRCh38 reference, including GL and KI supplementary contigs?

The first implementation was written around ordinary chromosome names such as `chr1`. GRCh38 also contains supplementary contigs whose identifiers include multiple periods and version suffixes, for example `.1` or `.2`.

The workflow needed to preserve each contig identifier exactly because even a small naming change breaks the connection among:

```text
reference FASTA
BAM contig
REDItools2 interval
coverage output
candidate coordinate
```

## Build

REDItools2 coverage generation was split by contig and later reconstructed. Temporary filenames encoded the original contig and interval information.

The relevant scripts include:

```text
pipeline/scripts/rna/generate_reditools_coverage_limited.sh
pipeline/scripts/rna/run_reditools_all_samples.sh
pipeline/scripts/rna/rebuild_reditools_file_list.py
```

The intended parsing rule was to remove only the compression suffix, then split encoded fields from the right.

## Test

### Failure around a supplementary contig

The workflow failed when processing a contig such as:

```text
chrGL000009.2
```

The original filename parser removed text based on the final period. That logic interpreted the contig version suffix as a file extension and changed the identifier.

A contig recorded in the FASTA as:

```text
chrGL000009.2
```

could incorrectly become:

```text
chrGL000009
```

The resulting identifier no longer matched the reference or BAM header.

### Corrected parser

The retained parser removes exactly the final three characters corresponding to `.gz`, then separates encoded fields from the right:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

This preserves the complete contig name, including `.1` or `.2`.

### Apparent stalls on GL/KI contigs

Coverage jobs sometimes appeared to stop when logs reached GL or KI contigs. These names were initially suspected to indicate an error.

Process and file-growth checks showed that the jobs could still be active:

```bash
ps -ef | grep samtools
find "$PROJECT/reditools" -type f -printf '%s %p\n' | sort -n | tail
```

The presence of a supplementary contig in a log is not by itself evidence of failure.

### Library warnings

`samtools` also emitted warnings involving Conda-provided ncurses or terminfo libraries, including references to:

```text
libncursesw.so.6
libtinfow.so.6
```

When valid depth output continued to be produced, the warning was treated as an environment warning rather than an automatic pipeline failure.

## Learn

### Lesson 1 — Biological identifiers are not filenames

Periods in contig names are part of the assembly identifier, not necessarily extension delimiters. Generic filename manipulation can silently alter biological coordinates.

### Lesson 2 — Full-reference support must be tested

Testing only canonical chromosomes would not have exposed this bug. The corrected workflow retains supplementary contigs and documents their expected appearance.

### Lesson 3 — Progress must be diagnosed from evidence

A long-running job is assessed using:

```text
active processes
output-file growth
log timestamps
exit status
output integrity
```

rather than the visual unfamiliarity of the current contig name.

### Lesson 4 — Warnings and failures are different classes

The revised workflow distinguishes:

- a warning with valid, growing output;
- a stalled process;
- a process that exited non-zero;
- a malformed or missing output file.

## Final safeguards

1. Preserve complete contig identifiers.
2. Remove only known compression suffixes.
3. Split encoded temporary filenames from the right.
4. Verify output growth before declaring a coverage job stalled.
5. Validate that output contigs exist in the reference index.

## Final role in the pipeline

This cycle ensured that the candidate-generation branch operates on the full GRCh38 reference rather than only canonical chromosomes, while preventing silent coordinate corruption during temporary-file handling.
