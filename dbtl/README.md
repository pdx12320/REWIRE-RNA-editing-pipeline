# Design–Build–Test–Learn record

This folder records how the RNA-editing screening workflow developed. The wiki presents the finished method and results; this folder keeps the iterative analysis, implementation changes, failures and decisions.

## Objective

Identify transcript-oriented C-to-U screening candidates that were called reproducibly in three treated libraries, not called in three controls under the original REDItools2 settings, consistent with transcript strand and absent from the harmonized 293T catalogue by exact allele.

The final set remains a screening set because the frozen legacy table does not contain independent control-site depth and base counts for control non-calls.

## Ten cycles

| Cycle | Question | Test | Decision |
|---|---|---|---|
| [1](cycle-1-rna-evidence.md) | Can RNA-seq provide reproducible editing evidence? | Process three treated and three control libraries independently | Retain replicate-aware RNA evidence |
| [2](cycle-2-public-wgs.md) | Can public WGS provide a genomic blacklist? | Align and filter three public WGS runs | Reject the route because mapping and variant yield were inadequate |
| [3](cycle-3-catalogue-harmonization.md) | Can `293T_CG` replace the failed WGS route? | Convert hg18 to GRCh38 and validate REF alleles | Adopt the catalogue with QC safeguards |
| [4](cycle-4-final-integration.md) | How should RNA and genomic evidence be combined? | Exact-allele comparison | Retain 3,333 catalogue-filtered screening candidates |
| [5](cycle-5-gatk-read-groups.md) | Can GATK preprocessing preserve sample identity? | Test read groups, BAM integrity and indexes | Validate or add read groups before duplicate handling |
| [6](cycle-6-reditools-environment.md) | Can REDItools2 run reproducibly with Python 2 and MPI? | Test dependencies, interpreter identity and worker allocation | Freeze a dedicated environment and pass the interpreter explicitly |
| [7](cycle-7-contigs-and-coverage.md) | Can supplementary GRCh38 contigs be handled safely? | Test GL/KI coverage jobs and filename parsing | Preserve complete contig identifiers |
| [8](cycle-8-control-depth-evidence.md) | Does a control non-call prove absence? | Separate caller output from direct candidate depth | Keep the control-depth boundary explicit |
| [9](cycle-9-legacy-compatibility.md) | How can the legacy table use the new catalogue? | Compare strict-schema requirements with available files | Add a narrow exact-allele compatibility route |
| [10](cycle-10-audit-and-reporting.md) | How should the result be reported without overclaiming? | Audit claims, counts and files | Preserve exclusions and freeze evidence-calibrated terminology |

## Supporting records

- [Failure log](failure-log.md)
- [Decision log](decision-log.md)
- [Reproducibility checklist](reproducibility-checklist.md)

## Final evidence funnel

```text
9,930 strand-consistent sites
→ 4,778 called in all three treated replicates
→ 3,349 treatment-specific candidates before catalogue comparison
→ 16 exact catalogue matches removed
→ 3,333 final screening candidates
```

## Repository relationship

```text
wiki/README.md       finished method and results
dbtl/                ten development cycles
pipeline/            executable implementation
results/              frozen summary
```
