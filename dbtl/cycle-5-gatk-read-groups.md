# Cycle 5 — Stabilise GATK RNA preprocessing

## Design

### Question

Could the alignment output be processed reproducibly by GATK without losing sample identity or splice-aware read structure?

The RNA-editing caller requires clean, coordinate-consistent BAM files. Before REDItools2, the alignments therefore had to pass through duplicate handling and splice-aware preprocessing.

### Initial plan

```text
STAR BAM
→ GATK MarkDuplicates
→ GATK SplitNCigarReads
→ sorted and indexed BAM
```

The intended output for each of the six libraries was one BAM that:

- retained sample identity;
- contained valid read-group metadata;
- was coordinate-sorted;
- had a matching BAM index;
- could be read by both GATK and samtools.

## Build

The preprocessing wrapper was implemented in:

```text
pipeline/scripts/rna/run_gatk_preprocessing.py
```

Typical invocation:

```bash
python3 pipeline/scripts/rna/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest pipeline/config/samples.tsv \
  --java-options=-Xmx16g
```

For each sample, the wrapper connects the expected STAR BAM, GATK output directories, metrics files and final indexes.

## Test

### Failure observed

GATK reported:

```text
SAMRecord.getReadGroup() is null
```

The BAM file could be opened, but the header did not provide the read-group information required by GATK.

### Diagnostic test

```bash
samtools view -H sample.bam | grep '^@RG'
```

The absence of an `@RG` line confirmed that the problem was metadata rather than damaged sequence alignments.

### Why this mattered

Without a valid read group, GATK cannot reliably associate reads with a sample. Continuing with an undefined sample identity would make downstream files difficult to audit and could cause errors in tools that rely on read-group fields.

## Learn

### Pipeline change

Read-group validation was moved before duplicate marking. When the STAR output lacks a usable read group, it must be added or replaced before `MarkDuplicates`.

The required metadata includes:

```text
ID    read-group identifier
SM    sample name
PL    sequencing platform
LB    library identifier
PU    platform unit, when available
```

### Validation rule retained

A BAM is not considered ready for GATK only because `samtools quickcheck` succeeds. It must also pass:

```bash
samtools view -H sample.bam | grep '^@RG'
samtools quickcheck sample.bam
```

After preprocessing, the final BAM and index are checked together:

```bash
samtools quickcheck sample.splitncigarreads.bam
test -s sample.splitncigarreads.bam.bai
```

### Broader lesson

File existence is weaker than workflow compatibility. The first implementation treated a readable BAM as sufficient; the revised workflow tests the metadata contract required by the next tool.

## Final role in the pipeline

This cycle established a stable handoff:

```text
STAR alignment
→ read-group validation
→ duplicate metrics
→ SplitNCigarReads
→ indexed BAM for REDItools2 and depth analysis
```

It also created a reusable rule for later steps: validate the properties required by the downstream tool, not only the file format itself.
