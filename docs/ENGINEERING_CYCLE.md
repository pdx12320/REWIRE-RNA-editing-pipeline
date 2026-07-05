# Engineering cycle and decision log

This document records how the final workflow was reached. The wiki keeps only the core story; the commands, failed tests and implementation decisions are preserved here.

## Final objective

Identify transcript-oriented C-to-U RNA-editing screening candidates that are:

1. called in all three treated RNA-seq libraries;
2. not called in the three controls under the original REDItools2 filter;
3. consistent with transcript strand;
4. absent from an assembly-harmonized 293T genomic-variant catalogue by exact `CHROM:POS:REF:ALT` matching.

The final output is a screening set, not a definitive off-target list, because the frozen legacy table does not contain independent depth and base counts for control non-calls.

---

# DBTL cycle 1 — Build the RNA evidence branch

## Design

Use three treated and three control RNA-seq libraries. Process each library independently so that replicate support remains visible.

## Build

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
STAR_INDEX=/data/ydx/igem/STAR_index
REDITOOLS=/data/ydx/igem/REDItools2
```

The implemented stages are:

```text
SRA Toolkit
→ STAR two-pass alignment
→ GATK MarkDuplicates
→ GATK SplitNCigarReads
→ REDItools2 per sample
→ union candidate coordinates
→ VEP transcript-strand annotation
→ treated/control comparison
```

Core scripts:

```text
pipeline/scripts/rna/download_sra_fastq.py
pipeline/scripts/rna/run_star_alignment.py
pipeline/scripts/rna/run_gatk_preprocessing.py
pipeline/scripts/rna/run_reditools_all_samples.sh
pipeline/scripts/rna/reditools_union_to_vcf.py
pipeline/scripts/rna/run_vep_annotation.py
pipeline/scripts/rna/filter_c_to_u_and_compare.py
```

## Test

The frozen run produced:

```text
strand-consistent site matrix             9,930
called in all three treated replicates    4,778
treatment-specific before catalogue       3,349
```

Transcript orientation was interpreted as:

```text
positive-strand transcript: genomic C>T
negative-strand transcript: genomic G>A
```

## Learn

A control non-call is not equivalent to zero editing. REDItools2 may omit a position because of its edited-read threshold. The strict route therefore includes direct candidate-site depth measurement in all six BAM files:

```bash
bash pipeline/scripts/rna/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  pipeline/config/samples.tsv
```

The frozen legacy output did not retain those independent control-depth values. We therefore report its retained sites as screening candidates.

---

# DBTL cycle 2 — Test public WGS as a genomic blacklist

## Design

The first genomic-variant strategy was to reconstruct a 293T blacklist from three public WGS runs.

## Build

Each run was evaluated separately, with a planned two-of-three exact-allele consensus when the runs represented different BioSamples.

## Test

Observed mapping rates were only 19.2–26.3%. Strict per-sample filtering retained 137–230 variants, and the two-of-three consensus contained only 118 variants.

## Learn

A blacklist with only 118 variants cannot support genome-wide exclusion. The public-SRA WGS route was removed from the final workflow rather than presented as successful validation.

Decision:

```text
abandon low-coverage WGS reconstruction
→ use the database-released 293T_CG variant catalogue
```

---

# DBTL cycle 3 — Harmonize the 293T_CG catalogue

## Design

Use the HEK293 Genome Project `293T_CG` VCF as external genomic evidence. The source catalogue uses NCBI build 36/hg18, while the RNA branch uses GRCh38, so coordinate conversion and reference validation are mandatory.

## Build

Environment:

```bash
conda env create -f pipeline/env/genomic_catalogue.yml
conda activate rewire_catalogue
```

Run:

```bash
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
CATALOGUE_IN=/data/ydx/igem/293T_CG.vcf
CHAIN=/data/ydx/igem/hg18ToHg38.over.chain.gz
CATALOGUE_OUT=/data/ydx/igem/293T_CG_GRCh38

bash pipeline/scripts/catalogue/process_293T_CG_to_GRCh38.sh \
  --input "$CATALOGUE_IN" \
  --reference "$REF" \
  --chain "$CHAIN" \
  --outdir "$CATALOGUE_OUT" \
  --threads 16
```

The script performs:

```text
content-based VCF format detection
→ PASS biallelic SNP selection
→ CrossMap hg18-to-GRCh38 liftover
→ remove unmapped records
→ GRCh38 REF validation
→ remove REF mismatches
→ normalize and coordinate-sort
→ bgzip and tabix index
→ QC summary and checksums
```

## Test

```text
source PASS biallelic SNPs       2,914,465
CrossMap-unmapped records            5,979
GRCh38 REF mismatches removed       22,761
final GRCh38 SNPs                2,885,725
```

## Learn

Several failures changed the implementation:

| Failure | Diagnosis | Change |
|---|---|---|
| `Exec format error` when reading `293T_CG.vcf` | the file content was gzip-compressed despite the extension | detect format from file bytes rather than filename |
| UCSC chain download timed out | network retrieval was unreliable on the server | accept a validated local `hg18ToHg38.over.chain.gz` file |
| tabix indexing failed with unsorted positions | CrossMap output was not coordinate-sorted | run `bcftools sort` before indexing |
| mapped records had incompatible REF alleles | liftover coordinates alone were insufficient | validate against the project GRCh38 FASTA with `bcftools norm -c x` |

These checks were retained in the final script because each one corresponds to an observed failure, not a hypothetical edge case.

---

# DBTL cycle 4 — Integrate the catalogue with the frozen RNA result

## Design

Match exact alleles rather than coordinates alone:

```text
CHROM : POS : REF : ALT
```

Keep excluded records in a separate table so that the filtering decision remains auditable.

## Build

The current strict route expects `candidate_depth/` and a site matrix containing `all_replicates_depth_pass`. The frozen legacy table lacked both. Instead of fabricating a depth flag, we added a compatibility script:

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
CATALOGUE=/data/ydx/igem/293T_CG_GRCh38_retry/293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz

python3 pipeline/scripts/catalogue/filter_existing_treatment_specific_by_293T.py \
  --treatment-specific "$PROJECT/final/CU5.17_EGFP_GC.treatment_specific.tsv.gz" \
  --catalogue-vcf "$CATALOGUE" \
  --output-dir "$PROJECT/final_with_293T_catalogue"
```

## Test

```text
treatment-specific before catalogue    3,349
exact catalogue overlaps                  16
retained screening candidates          3,333
```

Outputs:

```text
CU5.17_EGFP_GC.treatment_specific.before_catalogue.annotated.tsv.gz
CU5.17_EGFP_GC.catalogue_overlaps.tsv.gz
CU5.17_EGFP_GC.treatment_specific.tsv.gz
catalogue_integration_summary.tsv
```

## Learn

The final workflow should preserve the evidence boundary:

- the 16 overlapping alleles are plausible genomic variants and are excluded from the retained set;
- they are not negative editing examples;
- absence from the catalogue does not prove absence of a genomic variant in the exact experimental subline;
- the 3,333 retained records are screening candidates until control-site depth and base counts are measured directly.

---

# Final file map

```text
wiki/README.md
    concise iGEM-facing story

docs/ENGINEERING_CYCLE.md
    DBTL history, failed tests, decisions and commands

pipeline/README.md
    executable workflow

pipeline/scripts/rna/
    RNA processing and evidence integration

pipeline/scripts/catalogue/
    catalogue harmonization and compatibility filtering

pipeline/CATALOGUE_PROVENANCE.md
    source, assembly conversion and QC counts

pipeline/OUTPUTS.md
    expected files

pipeline/TROUBLESHOOTING.md
    observed failure modes and fixes

results/final_summary.tsv
    frozen evidence-funnel counts
```

# Reproducibility rules

1. Use the same GRCh38 FASTA across STAR, GATK, REDItools2, CrossMap validation and bcftools normalization.
2. Record software versions and file checksums.
3. Preserve supplementary-contig version suffixes.
4. Do not equate a missing control call with zero editing.
5. Match catalogue evidence by exact allele, not coordinate alone.
6. Keep excluded records for audit.
7. Do not describe the retained set as definitively SNV-free or fully depth-qualified.
