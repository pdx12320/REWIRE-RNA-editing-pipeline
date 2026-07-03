# Evidence filtering strategy

## Why a missing call is not enough

REDItools2 reports substitutions that satisfy its internal thresholds. When a site is absent from a control call table, two explanations remain possible: the control may truly lack the substitution, or the position may have insufficient sequencing coverage. We therefore measure candidate-site depth independently in every treated and control BAM.

## Strand rule

Transcript-level C-to-U editing is represented differently in genomic coordinates:

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

Only substitutions consistent with the VEP-annotated transcript orientation are retained. Coordinates with conflicting strand assignments are excluded.

## Default replicate rule

The conservative treatment-specific definition requires:

```text
called in 3 of 3 treated replicates
called in 0 of 3 control replicates
candidate-site depth ≥20 in all six samples
base quality ≥30 for the independent depth check
mapping quality ≥20 for the independent depth check
not present in the optional matched HEK293T WGS VCF
```

## Why the full matrix is retained

The pipeline does not export only a final binary label. For each candidate coordinate, the evidence matrix stores:

- treated and control replicate call counts;
- per-sample call status;
- per-sample candidate-site depth;
- REDItools2 quality-filtered depth when called;
- edited-read count;
- editing frequency;
- VEP transcript strand;
- optional WGS overlap status.

This design keeps the decision process auditable. Less stringent definitions, such as detection in two of three treated replicates, can be tested later without repeating alignment or REDItools2 calling.
