# Cycle 8 — Distinguish control non-calls from control-negative evidence

## Design

### Question

When a candidate is absent from a control REDItools2 table, does that mean the control truly lacks editing at that site?

The first screening rule used call status:

```text
called in all three treated replicates
AND not called in any control replicate
```

This is useful for candidate reduction, but REDItools2 was run with an edited-read threshold. A site can therefore be missing from the control call table for two very different reasons:

1. the control had sufficient coverage and little or no alternate signal;
2. the control had insufficient depth or did not reach the reporting threshold.

The pipeline needed an independent coverage layer to separate these cases.

## Build

### Candidate union

All reported alleles were first merged into a union BED:

```text
$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed
```

### Direct depth measurement

Candidate depth was measured from each final SplitNCigarReads BAM rather than inferred from whether REDItools2 emitted a call.

```bash
bash pipeline/scripts/rna/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  pipeline/config/samples.tsv
```

Expected outputs:

```text
$PROJECT/candidate_depth/CU517_GC_T1.candidate_depth.tsv.gz
$PROJECT/candidate_depth/CU517_GC_T2.candidate_depth.tsv.gz
$PROJECT/candidate_depth/CU517_GC_T3.candidate_depth.tsv.gz
$PROJECT/candidate_depth/CU517_GC_C1.candidate_depth.tsv.gz
$PROJECT/candidate_depth/CU517_GC_C2.candidate_depth.tsv.gz
$PROJECT/candidate_depth/CU517_GC_C3.candidate_depth.tsv.gz
```

Each record contains:

```text
chromosome
1-based genomic position
depth
```

### Strict integration rule

The current integration route can require:

```text
minimum depth >=20 in all six libraries
```

through:

```text
all_replicates_depth_pass
```

This creates a stricter definition:

```text
3/3 treated calls
0/3 control calls
adequate depth in all six libraries
strand-consistent C-to-U
```

## Test

### Frozen result behavior

The earlier completed workflow produced:

```text
site_matrix = 9,930
treated_consensus = 4,778
treatment_specific = 3,349
```

However, during final catalogue integration, the available legacy site matrix did not contain:

```text
all_replicates_depth_pass
```

and the original `candidate_depth/` directory was not available at the expected location.

### Failed assumption avoided

It would have been technically easy to:

- rename another field;
- fill the missing depth flag with `1`;
- treat every control non-call as zero editing.

All three choices would create evidence that had not actually been measured. The pipeline therefore stopped rather than silently continuing.

## Learn

### Lesson 1 — Caller output and BAM evidence are different layers

A call table answers:

> Did this site pass the caller's reporting criteria?

A depth table answers:

> Was this position observable in this sample?

These are not interchangeable.

### Lesson 2 — Missing evidence must remain visible

The frozen legacy output remains useful for screening, but its interpretation is limited. The final 3,333 retained records are therefore not described as fully depth-qualified control-negative events.

### Lesson 3 — The strict route remains available

For a future rerun, the complete workflow should regenerate candidate-depth tables and apply the all-replicate depth gate before defining the final high-confidence set.

### Lesson 4 — Thresholds should be recorded separately

The project contains several thresholds with different meanings:

| Threshold | Meaning |
|---|---|
| REDItools2 `-me 20` | minimum edited-read support for a reported call |
| mapping quality | read-alignment confidence |
| base quality | confidence in the observed nucleotide |
| candidate depth ≥20 | independent observability at the candidate coordinate |

Combining these into one unnamed “quality filter” would hide the logic of the result.

## Final reporting rule

For the frozen output, use:

> called in all three treated replicates and not called in the controls under the original REDItools2 settings

Do not replace this with:

> absent from all controls

unless direct depth and base-count evidence has been generated and reviewed.

## Final role in the pipeline

This cycle established the evidence hierarchy used throughout the repository:

```text
call status
< independent candidate depth
< direct base counts and matched DNA
```

It also explains why the final result is presented as a screening set with an explicit control-evidence limitation.
