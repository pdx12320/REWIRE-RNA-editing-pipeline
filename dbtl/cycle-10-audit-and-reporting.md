# Cycle 10 — Make the final result auditable and evidence-calibrated

## Design

### Question

How should the final dry-lab result be reported so that another team can understand what was measured, reproduce the filtering, and avoid overstating the evidence?

A technically correct pipeline can still produce a misleading project page when:

- excluded sites disappear without explanation;
- several distinct thresholds are compressed into “quality filtering”;
- an external catalogue is described as matched WGS;
- control non-calls are described as confirmed zero editing;
- candidate sites are called validated off-targets before orthogonal testing.

The final reporting layer therefore became part of the engineering design rather than an afterthought.

## Build

### Output separation

The final catalogue integration writes four outputs with different roles:

```text
1. annotated pre-catalogue table
2. catalogue-overlap exclusion table
3. retained screening table
4. count summary
```

This structure preserves both the retained evidence and the filtering history.

### Frozen summary

```text
strand-consistent site matrix             9,930
called in all three treated replicates    4,778
treatment-specific before catalogue       3,349
exact 293T catalogue overlaps                16
final catalogue-filtered candidates       3,333
```

The compact machine-readable version is stored in:

```text
results/final_summary.tsv
```

### Repository separation

The public-facing and technical records were separated by purpose:

```text
wiki/README.md
    finished method, results, contribution and limitations

dbtl/
    iterative design history and lessons

pipeline/
    executable implementation and commands

results/
    frozen result summary
```

The wiki does not need every server failure. The DBTL folder does not need to repeat the polished project narrative. Both remain linked.

### Terminology ledger

The following terms were fixed across the repository:

| Term | Meaning |
|---|---|
| `treated consensus` | called in all three treated replicates under the implemented call rules |
| `treatment-specific before catalogue` | treated consensus and not called in controls under the original settings |
| `genomic_catalogue_overlap` | exact `CHROM:POS:REF:ALT` present in the external 293T catalogue |
| `final screening candidates` | treatment-specific candidates after exact catalogue exclusions |

Preferred final phrase:

> catalogue-filtered treatment-associated C-to-U screening candidates

## Test

### Claim audit

Each candidate claim was compared with the available evidence.

#### Claim: “The site is an RNA-editing off-target”

Evidence available:

```text
RNA-seq mismatch
replicate call support
control-call comparison
transcript orientation
external catalogue comparison
```

Evidence missing:

```text
matched DNA from the exact experimental cells
independent control-site base counts in the frozen table
orthogonal amplicon validation
```

Decision: do not use “confirmed off-target.”

#### Claim: “The site is absent from controls”

Evidence available: no REDItools2 control call.

Evidence missing: independent depth and base counts for the frozen output.

Decision: write “not called in controls under the original REDItools2 settings.”

#### Claim: “The site is not a genomic variant”

Evidence available: absent from one harmonized external 293T catalogue.

Evidence missing: matched WGS for the experimental subline.

Decision: write “did not overlap the selected 293T catalogue.”

### Count audit

The final outputs were checked for internal consistency:

```text
retained rows + excluded rows = annotated pre-catalogue rows
3,333 + 16 = 3,349
```

The summary table and narrative numbers were aligned to the same frozen values.

### File audit

Large raw data were kept outside GitHub, while the repository retained:

- scripts;
- environment definitions;
- sample manifest;
- catalogue provenance;
- output contracts;
- frozen counts;
- decision and failure logs.

This balance keeps the repository usable without pretending that raw sequencing data belong in source control.

## Learn

### Lesson 1 — Scientific wording is part of model validation

Overstated language can make a result appear stronger than the data. The final terminology was therefore treated as a formal boundary of the pipeline.

### Lesson 2 — Exclusions are results

The 16 catalogue-overlapping alleles are not discarded noise. They demonstrate that the catalogue comparison changed the candidate set and provide examples for auditing exact-allele behavior.

### Lesson 3 — Reproducibility requires both code and rationale

Code shows how a result was produced. DBTL records explain why the final route was selected over alternatives. Both are required for reuse.

### Lesson 4 — The public wiki should follow reader order

The final wiki presents:

```text
biological problem
→ input design
→ core workflow
→ result funnel
→ interpretation
→ contribution
→ limitations
```

Detailed installation and failure history remain linked but do not interrupt the main dry-lab story.

## Final contribution of this cycle

This cycle produced an evidence-calibrated deliverable rather than only a list of coordinates:

- retained candidates;
- explicit exclusions;
- frozen counts;
- reproducible code;
- catalogue provenance;
- ten DBTL cycles;
- stated limitations and future validation requirements.

The final output is therefore suitable for iGEM reporting and for transfer to another analyst without hiding the uncertainty that remains.
