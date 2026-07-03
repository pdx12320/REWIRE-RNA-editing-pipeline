# REWIRE RNA-editing pipeline

This repository contains the paper-style workflow currently used to analyze the CU5.17_EGFP_GC RNA-seq dataset.

## Workflow

```text
SRA paired-end RNA-seq
→ STAR two-pass alignment to GRCh38
→ GATK read-group validation
→ MarkDuplicates
→ SplitNCigarReads
→ REDItools2 parallel calling
→ VEP transcript-strand annotation
→ transcript-oriented C-to-U filtering
→ treated/control replicate comparison
→ optional HEK293T WGS variant removal
```

## Dataset

| Group | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

The machine-readable manifest is stored in `config/samples.tsv`.

## Repository layout

```text
config/   fixed sample manifest
docs/     installation, runbook, workflow, outputs, troubleshooting, limitations
scripts/  REDItools2 coverage, union-VCF, strand filtering, replicate comparison
wiki/     English text for direct adaptation to an iGEM wiki
results/  placeholders for final tables and figures
```

## Main scripts

- `generate_reditools_coverage_limited.sh` builds per-sample coverage maps without launching one disk-intensive job for every contig simultaneously.
- `reditools_union_to_vcf.py` combines substitutions from the six REDItools2 tables into a union VCF and candidate-position BED.
- `filter_c_to_u_and_compare.py` applies transcript-strand interpretation, all-replicate depth checks, treated-control comparison, and optional WGS filtering.

Run shell files with `bash` and Python files with `python3`.

## Methodological notes

REDItools2 is run in strict mode with a minimum of 20 edited reads. This is an edited-read threshold, not simply a total-depth cutoff. Candidate-site depth is measured separately in all six replicates using base quality 30 and mapping quality 20.

Transcript-level C-to-U editing is represented as genomic C-to-T on positive-strand transcripts and genomic G-to-A on negative-strand transcripts. Conflicting transcript-strand annotations are excluded.

Without a matched HEK293T WGS variant file, the final table remains an RNA-derived candidate list rather than a definitive editing-only set.

## Results status

The analysis is still running. No numerical result counts from the older lightweight CU5.15 analysis are retained on the main branch. Final CU5.17 tables and figures will be added only after all six samples complete the same workflow and pass integrity checks.

## Documentation

- `docs/INSTALLATION.md`
- `docs/RUNBOOK.md`
- `docs/PIPELINE.md`
- `docs/OUTPUTS.md`
- `docs/TROUBLESHOOTING.md`
- `docs/LIMITATIONS.md`
- `wiki/Dry_Lab.md`
- `wiki/Methods.md`
- `wiki/RESULTS_TEMPLATE.md`
