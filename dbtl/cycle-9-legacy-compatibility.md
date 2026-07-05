# Cycle 9 — Add a compatibility route for the completed legacy output

## Design

### Question

How could the completed 3,349-row treatment-specific table be compared with the new 293T catalogue when the current strict integration inputs were incomplete?

By the time the catalogue branch was ready, the RNA branch had already been completed with an earlier helper implementation. The newer integration code expected:

```text
candidate_depth/ directory
site matrix with all_replicates_depth_pass
current filter_calls.py and filter_utils.py schema
```

The available legacy result contained the final treatment-specific records but not the full current schema.

The goal was to preserve the valid operation still supported by the data—exact catalogue comparison—without reconstructing or inventing missing depth evidence.

## Build

### Initial failure

The first final-integration wrapper stopped with:

```text
Required directory is missing: .../candidate_depth
```

After attempting to use the existing site matrix directly, the integration also reported:

```text
site matrix is missing required columns: all_replicates_depth_pass
```

These failures showed that the strict current route and the frozen legacy output represented different data contracts.

### Alternatives considered

#### Option A — Fill the missing field

Rejected because no evidence supported assigning every row a passing depth status.

#### Option B — Copy only one updated filter script

Rejected because `filter_c_to_u_and_compare.py`, `filter_calls.py` and `filter_utils.py` form one coordinated implementation. Mixing old and new files risks schema and terminology mismatches.

#### Option C — Rerun the entire RNA pipeline

Possible in principle, but unnecessary for the narrow catalogue-comparison task and dependent on reconstructing all original intermediate files.

#### Option D — Build a narrow compatibility filter

Selected because the existing 3,349-row treatment-specific table already represented the completed RNA screening output. The remaining supported operation was exact comparison against the external catalogue.

### Compatibility script

```text
pipeline/scripts/catalogue/filter_existing_treatment_specific_by_293T.py
```

Inputs:

```text
completed treatment-specific TSV or TSV.GZ
GRCh38 catalogue VCF or VCF.GZ
output directory
```

The script:

1. detects chromosome, position, REF and ALT columns;
2. normalizes `chr` prefixes without changing contig identity;
3. loads PASS catalogue alleles;
4. compares exact `CHROM:POS:REF:ALT` keys;
5. adds `genomic_catalogue_overlap`;
6. writes annotated, excluded and retained tables;
7. writes a count summary.

Frozen command:

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
INPUT="$PROJECT/final/CU5.17_EGFP_GC.treatment_specific.tsv.gz"
CATALOGUE=/data/ydx/igem/293T_CG_GRCh38_retry/293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz
OUT="$PROJECT/final_with_293T_catalogue"

python3 pipeline/scripts/catalogue/filter_existing_treatment_specific_by_293T.py \
  --treatment-specific "$INPUT" \
  --catalogue-vcf "$CATALOGUE" \
  --output-dir "$OUT"
```

## Test

The script loaded:

```text
2,885,725 catalogue exact alleles
```

and produced:

```text
treatment_specific_before_catalogue = 3,349
catalogue_overlap                    = 16
final_treatment_specific             = 3,333
```

Output contract:

```text
CU5.17_EGFP_GC.treatment_specific.before_catalogue.annotated.tsv.gz
CU5.17_EGFP_GC.catalogue_overlaps.tsv.gz
CU5.17_EGFP_GC.treatment_specific.tsv.gz
catalogue_integration_summary.tsv
```

### Arithmetic validation

```text
3,349 - 16 = 3,333
```

This simple subtraction was included as a consistency check between annotated, excluded and retained outputs.

## Learn

### Lesson 1 — Version compatibility should be explicit

A pipeline update can change required columns and intermediate files. Compatibility should be handled by a named route rather than hidden inside ad hoc file edits.

### Lesson 2 — Narrow scripts can preserve validity

The compatibility script does not pretend to recreate the strict full workflow. It performs one defined operation supported by the available table.

### Lesson 3 — Old results require frozen semantics

The meaning of the 3,349 input rows remains tied to the original filtering implementation. Catalogue filtering adds genomic-overlap evidence but does not strengthen the missing control-depth evidence.

### Lesson 4 — Schema checks prevent silent errors

The script explicitly searches for required allele columns and stops when they are absent. It does not assume a fixed column position.

## Final role in the pipeline

The repository now supports two routes:

```text
Route A: strict current integration
    complete REDItools + VEP + candidate-depth inputs

Route B: legacy compatibility integration
    completed treatment-specific table + catalogue
```

This distinction allows the frozen result to remain reproducible while preserving a stricter route for future reruns.
