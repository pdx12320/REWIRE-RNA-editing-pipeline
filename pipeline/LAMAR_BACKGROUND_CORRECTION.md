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
or independent Lamar training set. It is a subset of the broad matrix and must
not be used as an independent test set when broad-matrix sites are used for
training.

## Intended workflow versus the actual frozen run

The repository's intended preprocessing route runs STAR and duplicate handling
consistently for every sample. The completed 2026-07-15 audit used the files that
were actually available on the server:

| Samples | Frozen BAM input | Duplicate history |
|---|---|---|
| T1 | Coordinate-sorted Picard MarkDuplicates BAM | Duplicates marked, not removed. |
| T2, T3, C1, C2, C3 | Original STAR coordinate-sorted BAMs | Duplicate marking had not been applied. |

The pileup excluded reads carrying the duplicate flag consistently in all six
inputs. Nevertheless, the five STAR BAMs did not contain Picard duplicate marks,
so the preprocessing histories were not identical. This limitation is frozen in
`bam_qc.tsv`, `repository_audit.md` and the analysis summary.

`training_eligible=1` requires the frozen sequence/orientation checks and at
least 2 of 3 sufficiently covered replicates in both treated and control groups.
`label_confidence=high` is stricter: all 3+3 replicates must be covered and the
treated/control MAD consistency thresholds must pass. These definitions are not
changed by the downstream handoff builder.

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

## Software-version boundary

The frozen production manifest records Python 3.14.6 (free-threading,
conda-forge), pysam 0.24.0, samtools 1.23.1 and SciPy 1.18.0. Those exact values
are preserved in `pipeline/env/lamar_labels_production_versions.txt`.

Repository code and synthetic handoff tests support Python 3.11 and 3.12 in CI.
The maintained BAM environment pins the currently evidenced pysam 0.24.x and
samtools 1.23.x families; end-to-end production evidence is specifically at
pysam 0.24.0 and samtools 1.23.1. Python 3.11/3.12 support is not a claim that
those interpreters exactly reproduce the Python 3.14 production runtime.

SciPy is optional in the audited statistics code. When available,
`scipy.stats.fisher_exact` is used; otherwise the repository uses its tested
stdlib two-sided exact-test fallback. The optional scalar baseline uses
scikit-learn 1.7.x from `pipeline/env/lamar_scalar_baseline.yml`.

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
the repository's overlap-cluster split, which merges position±50 genomic
intervals and also groups identical sequences. The supported chromosome-held-out
strategy measures a stronger genomic distribution shift.

The validated derived handoff contains 9,428 all-eligible rows, including 1,564
valid zero corrected labels; 8,540 rows are high confidence and are recommended
for the primary analysis. The 7,351 high-confidence rows without the elevated
control-background flag form a stricter sensitivity analysis.

```bash
python pipeline/scripts/rna/prepare_lamar_finetuning_handoff.py \
  --labels /path/to/background_corrected_labels.tsv.gz \
  --metadata /path/to/lamar_ready_metadata.tsv.gz \
  --output-dir /path/to/CU5.17_lamar_finetuning_handoff \
  --seed 20260715 \
  --split-strategy overlap_cluster
```

This packaging step validates and partitions frozen derived tables only. It does
not rerun BAM pileups or alter scientific counts and thresholds.
