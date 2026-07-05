# Design–Build–Test–Learn record

This folder documents how the final RNA-editing screening workflow was developed. It is intentionally separated from the iGEM-facing wiki page: the wiki presents the finished method and results, whereas this folder records failed approaches, implementation changes, diagnostic evidence and the reasoning behind each decision.

## Project objective

The dry-lab task was to identify transcript-oriented C-to-U RNA-editing screening candidates that were:

1. detected reproducibly in all three treated RNA-seq libraries;
2. not called in the three controls under the original REDItools2 settings;
3. consistent with transcript strand;
4. absent from an assembly-harmonized 293T genomic-variant catalogue by exact `CHROM:POS:REF:ALT` matching.

The final output is a screening set rather than a definitive off-target list. The frozen legacy result table does not retain independent control-site depth and base counts for control non-calls.

## Cycle map

| Cycle | Initial question | Main test | Decision |
|---|---|---|---|
| [Cycle 1](cycle-1-rna-evidence.md) | Can RNA-seq provide reproducible transcript-oriented editing evidence? | Three treated and three control libraries processed independently | Retain replicate-aware RNA evidence and state the control-depth boundary |
| [Cycle 2](cycle-2-public-wgs.md) | Can public WGS provide a genomic blacklist? | Three public WGS runs aligned and filtered | Reject the route because mapping and variant yield were inadequate |
| [Cycle 3](cycle-3-catalogue-harmonization.md) | Can the HEK293 `293T_CG` catalogue replace the failed WGS route? | hg18-to-GRCh38 conversion, REF validation and QC | Adopt the catalogue after adding format, sorting and REF safeguards |
| [Cycle 4](cycle-4-final-integration.md) | How should RNA evidence and genomic evidence be combined? | Exact-allele comparison and audit-table generation | Retain 3,333 catalogue-filtered screening candidates |

## Supporting records

- [Failure log](failure-log.md): observed errors, diagnoses, fixes and retained safeguards.
- [Decision log](decision-log.md): major methodological decisions and the evidence supporting them.
- [Reproducibility checklist](reproducibility-checklist.md): files, versions, reference rules and checks required to repeat the analysis.

## Final evidence funnel

```text
9,930 strand-consistent sites
→ 4,778 called in all three treated replicates
→ 3,349 treatment-specific candidates before catalogue comparison
→ 16 exact matches to the 293T catalogue removed
→ 3,333 final catalogue-filtered screening candidates
```

## Relationship to the rest of the repository

```text
wiki/README.md
    concise finished method, results, contribution and limitations

dbtl/
    design history, failed tests, changes and lessons

pipeline/README.md
    executable workflow and commands

pipeline/scripts/rna/
    RNA processing and evidence integration

pipeline/scripts/catalogue/
    catalogue harmonization and exact-allele filtering

results/final_summary.tsv
    frozen evidence-funnel counts
```
