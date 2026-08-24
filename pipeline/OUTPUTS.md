# Expected outputs

## RNA-seq branch

```text
$PROJECT/
├── fastq/
├── bam/
│   ├── star/
│   ├── markduplicates/
│   └── splitncigarreads/
├── metrics/
├── reditools/
│   ├── coverage/
│   ├── tmp/
│   └── tables/
├── vcf/
├── vep/
├── candidate_depth/          optional strict-depth route
├── final/
├── lamar_training/           Model 1 to Model 2 label route
└── logs/
```

### Per-sample REDItools2 tables

```text
reditools/tables/CU517_GC_T1.txt.gz
reditools/tables/CU517_GC_T2.txt.gz
reditools/tables/CU517_GC_T3.txt.gz
reditools/tables/CU517_GC_C1.txt.gz
reditools/tables/CU517_GC_C2.txt.gz
reditools/tables/CU517_GC_C3.txt.gz
```

Each file should have the expected REDItools2 header and a tabix index when required by downstream tools.

## Strict integration outputs

```text
final/CU5.17_EGFP_GC.site_matrix.tsv.gz
final/CU5.17_EGFP_GC.treated_consensus.tsv.gz
final/CU5.17_EGFP_GC.treatment_specific.tsv.gz
```

The strict route retains sample-level call status, independent candidate depth, REDItools2 coverage, alternate-read count, editing fraction and `genomic_catalogue_overlap`.

## LAMAR training-label outputs

```text
lamar_training/
├── base_counts/
│   ├── CU517_GC_T1.candidate_base_counts.tsv.gz
│   ├── CU517_GC_T2.candidate_base_counts.tsv.gz
│   ├── CU517_GC_T3.candidate_base_counts.tsv.gz
│   ├── CU517_GC_C1.candidate_base_counts.tsv.gz
│   ├── CU517_GC_C2.candidate_base_counts.tsv.gz
│   └── CU517_GC_C3.candidate_base_counts.tsv.gz
└── CU5.17_EGFP_GC.lamar_training_labels.tsv.gz
```

Each per-sample table contains quality-filtered A/C/G/T counts, forward/reverse counts, ref/alt counts, allele depth, editing rate and exclusion diagnostics. The combined table retains all six sample measurements and adds transcript-oriented sequence windows, treated/control medians, corrected editing efficiency, replicate MAD, pooled Fisher/FDR, label class, confidence and training eligibility.

Prefer the broad strand-consistent site matrix as input. The final 3,333-site screening table is already selected for treated presence and control non-calls and should not be the sole Model 2 training set.

## Audited background-correction outputs

The production audit route writes immutable timestamped directories under
`lamar_background_corrected/` and updates `latest` only after all validations and
checksums pass. See [`LAMAR_BACKGROUND_CORRECTION.md`](LAMAR_BACKGROUND_CORRECTION.md)
for the complete file manifest, recovery behavior and frozen QC counts.

The frozen audit used a Picard MarkDuplicates BAM for T1 and original STAR
coordinate-sorted BAMs for T2/T3/C1/C2/C3. Duplicate-flagged reads were excluded
consistently, but the preprocessing histories were not identical.

## LAMAR fine-tuning handoff outputs

```text
CU5.17_lamar_finetuning_handoff/
├── CU5.17_lamar_all_eligible.tsv.gz
├── CU5.17_lamar_high_confidence.tsv.gz
├── CU5.17_lamar_high_confidence_low_control.tsv.gz
├── CU5.17_lamar_excluded.tsv.gz
├── CU5.17_lamar_splits.tsv.gz
├── data_dictionary.tsv
├── split_qc.json
├── handoff_manifest.json
├── checksums.sha256
└── README.md
```

Use the high-confidence table as the primary dataset. Use all eligible rows and
the stricter low-control-background table as sensitivity analyses. The split
table covers every eligible row and guarantees that no allele key, overlapping
101-nt genomic cluster or identical sequence crosses splits.

`data/processed/` contains only the compact public model-facing subset of these
files. BAM, FASTQ, FASTA, full pileups and server-specific symlinks remain
outside Git.

## Legacy-table compatibility outputs

```text
final_with_293T_catalogue/
├── CU5.17_EGFP_GC.treatment_specific.before_catalogue.annotated.tsv.gz
├── CU5.17_EGFP_GC.treatment_specific.tsv.gz
├── CU5.17_EGFP_GC.catalogue_overlaps.tsv.gz
└── catalogue_integration_summary.tsv
```

Frozen result:

```text
treatment-specific before catalogue    3,349
exact catalogue overlaps                  16
retained screening candidates          3,333
```

The 16 excluded alleles are preserved for audit. They should not be treated as negative editing examples.

## 293T catalogue branch

```text
$CATALOGUE_OUT/
├── 293T_CG.hg18.PASS.biallelic.SNV.vcf.gz
├── 293T_CG.hg18.PASS.biallelic.SNV.vcf.gz.tbi
├── 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
├── 293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz.tbi
├── logs/
├── qc/
└── tmp/
```

The frozen integration used:

```text
293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
```

Catalogue QC files:

```text
qc/input_header.txt
qc/summary.txt
qc/293T_CG.GRCh38.stats.txt
qc/variants_per_contig.txt
qc/first_20_variants.tsv
qc/checksums.sha256
```

The frozen conversion contains 2,885,725 final GRCh38 PASS biallelic SNPs.
