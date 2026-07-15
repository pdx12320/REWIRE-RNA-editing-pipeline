# Cycle 12 — Execute and audit six-sample background correction

## Design

### Question

Does the Model 1 to Model 2 interface remain valid when run on every broad
candidate with independently measured treated and control counts?

### Risks identified

- a control non-call could be silently converted to zero;
- chromosome or coordinate mismatches could return zero coverage;
- T1 preprocessing differed from the other five BAMs;
- an interrupted run could lose completed in-memory pileups;
- a successful script exit could be mistaken for biological validation.

## Build

The audited runner validates FASTA, BAM/index and candidate compatibility before
counting. It records the T1 selection, uses consistent read filters, preserves
coverage status and raw counts, requires sufficient coverage in both groups,
extracts transcript-oriented sequence context, and writes a timestamped manifest,
reports, commands, versions and checksums.

A strict `--reuse-pileup` loader was added after a power interruption occurred
after all six BAMs had been counted. Recovery is allowed only for a complete,
unique and internally consistent 59,580-row table.

## Test

- 7/7 unit tests passed.
- FASTA, BAM/index, GRCh38 1-based, chromosome, allele, duplicate and strand
  checks passed.
- Broad and final row counts were 9,930 and 3,333.
- Corrected efficiencies were never negative; raw differences retained negative
  values.
- No insufficient-control site received a corrected label or training eligibility.
- Twenty deterministic sites were recounted in all six BAMs with independent
  `samtools mpileup`: 120/120 comparisons passed.
- All output SHA-256 checksums passed before `latest` was updated.

## Learn

### Lesson 1 — descriptive statistics are not labels

A median computed from one covered control replicate may remain as an auditable
description, but the raw difference and corrected label must stay missing until
the required number of control replicates is covered.

### Lesson 2 — checkpoint recovery needs validation

Reusing a pileup table is safe only after validating its schema, samples, alleles,
uniqueness, row count and depth/base-count identities. A partial file must abort.

### Lesson 3 — preprocessing exceptions must remain visible

No original STAR T1 BAM was found. The valid coordinate-sorted Picard
MarkDuplicates BAM was used, with duplicates marked rather than removed. This
does not make the six BAM preprocessing histories identical.

### Lesson 4 — model-ready is format-specific

The output is ready for center-site scalar regression. Native Lamar token-level
training still requires a documented adapter, the true PUF target sequence and a
mask/weighting policy.

## Decision

Use the broad 9,930-site matrix as the primary audited label universe, retain the
3,333 post-catalogue table as a separate analysis, and expose only 9,428 broad
training-eligible rows to downstream scalar fine-tuning by default.
