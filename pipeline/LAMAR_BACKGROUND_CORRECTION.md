# Audited Lamar background-correction run

This route is the production/audit complement to the composable scripts in
[`LAMAR_TRAINING_LABELS.md`](LAMAR_TRAINING_LABELS.md). It was exercised on the
complete CU5.17 dataset and publishes a timestamped run only after all input,
count, label, sequence and direct-recount checks pass.

## What it adds

- verifies or creates the FASTA `.fai` and records all chromosome lengths;
- audits candidate schemas, GRCh38 1-based coordinates, reference alleles,
  chromosome naming, duplicate alleles and strand-specific C>T/G>A rules;
- searches recursively for the original STAR T1 BAM and records the actual T1
  choice instead of claiming identical preprocessing;
- independently counts A/C/G/T in T1, T2, T3, C1, C2 and C3 with MAPQ 30,
  base quality 20 and consistent flag filters;
- distinguishes adequate zero-alt coverage, low coverage, missing data and
  technical failure;
- requires sufficient treated and control coverage before assigning a corrected
  label or `training_eligible=1`;
- retains raw rates, replicate counts, medians, means, MADs, ranges, pooled counts,
  Fisher screening p-values and BH FDR;
- extracts a transcript-oriented 101-nt window and verifies that the center is C;
- validates 20 deterministic sites in all six BAMs with independent
  `samtools mpileup` calls (120 comparisons);
- writes checksums and atomically updates `latest` only after every validation.

Fisher exact testing is screening-only: reads are not independent biological
replicates. A successful run is computational QC, not experimental validation.

## Inputs selected by the audited CU5.17 route

The script expects the established project layout and selects:

```text
Broad matrix: final/CU5.17_EGFP_GC.site_matrix.tsv.gz
Final table:  final_with_293T_catalogue/CU5.17_EGFP_GC.treatment_specific.tsv.gz
```

The broad 9,930-site matrix is the primary label universe. The post-catalogue
3,333-site table is processed separately and must not be assumed to be a balanced
or independent Lamar training set.

## Run

```bash
conda env create -f pipeline/env/lamar_labels.yml
conda activate rewire_lamar_labels

PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa

nohup bash pipeline/scripts/rna/run_audited_lamar_background_correction.sh \
  "$PROJECT" "$REF" "$(command -v samtools)" "$(command -v python)" \
  > "$PROJECT/lamar_background_corrected/launcher.log" 2>&1 &
```

Do not start a second process while one is active. Monitor with:

```bash
tail -f "$PROJECT/lamar_background_corrected/launcher.log"
```

## Safe recovery from an interrupted post-pileup run

If all `9,930 × 6 = 59,580` broad pileup rows were written before an interruption,
the table can be re-used without rescanning BAMs:

```bash
bash pipeline/scripts/rna/run_audited_lamar_background_correction.sh \
  "$PROJECT" "$REF" "$(command -v samtools)" "$(command -v python)" \
  --reuse-pileup /path/to/incomplete_run/six_sample_pileup_counts.tsv.gz
```

The recovery loader aborts on missing columns, unexpected samples or alleles,
duplicate sample/allele rows, incorrect row count, or depth/base-count mismatch.
The reused table is recorded in `input_manifest.tsv`.

## Output

```text
$PROJECT/lamar_background_corrected/
├── run_<UTC timestamp>/
│   ├── input_manifest.tsv
│   ├── software_versions.txt
│   ├── commands.log
│   ├── repository_audit.md
│   ├── bam_qc.tsv
│   ├── reference_qc.txt
│   ├── candidate_table_qc.tsv
│   ├── selected_candidate_source.md
│   ├── six_sample_pileup_counts.tsv.gz
│   ├── background_corrected_labels.tsv.gz
│   ├── lamar_ready_metadata.tsv.gz
│   ├── excluded_sites.tsv.gz
│   ├── direct_pileup_validation.tsv
│   ├── qc_summary.tsv
│   ├── analysis_summary.md
│   ├── run.log
│   ├── checksums.sha256
│   └── final_candidate/
└── latest -> run_<UTC timestamp>
```

## Frozen validated run

The 2026-07-15 production run passed all checks:

| Metric | Broad | Final |
|---|---:|---:|
| Total sites | 9,930 | 3,333 |
| Training eligible | 9,428 | 3,328 |
| Sufficient coverage in all six | 9,091 | 3,327 |
| Positive raw differences | 7,864 | 3,327 |
| Zero after correction | 1,564 | 1 |
| Low treated coverage | 402 | 1 |
| Low control coverage | 356 | 5 |
| Elevated control background | 2,325 | 8 |
| BH FDR ≤ 0.05 | 7,412 | 3,317 |

All checksums passed. Independent recounting matched 120/120 sample-site
comparisons. The exact frozen counts are in
[`../results/lamar_background_corrected_qc_summary.tsv`](../results/lamar_background_corrected_qc_summary.tsv).

## Lamar fine-tuning boundary

`lamar_ready_metadata.tsv.gz` supplies one audited scalar center-site label per
101-nt sequence. The historical Lamar exporter uses
`seq,puf_target_seq,label_edit_rate,label_total_count` with 101-element token
arrays. A token-level trainer therefore needs an explicit center-mask adapter,
the experimentally correct PUF target sequence, and a documented weighting rule.
The PUF target must not be inferred from the sample name.

Avoid random row-level splits because nearby 101-nt windows can overlap. Prefer
chromosome-, gene-, transcript- or genomic-cluster-grouped evaluation.
