# Lamar handoff from REWIRE Model 1

## Frozen input set

The current handoff contains **3,333 unique GRCh38 candidate alleles** after
removing 16 exact matches to the assembly-harmonized 293T genomic catalogue.

All retained records satisfy the existing Model 1 treatment/control criteria
and have `genomic_catalogue_overlap = 0`.

- Positive-strand candidates are genomic C→T.
- Negative-strand candidates are genomic G→A and must be reverse-complemented
  so that Lamar receives transcript-level C-to-U context.
- `median_treated_edit_rate` is the recommended **provisional** continuous
  label because it is robust across the three treated replicates.
- `sd_treated_edit_rate` and `range_treated_edit_rate` quantify replicate
  disagreement.

## Generate the handoff table

Without sequence extraction:

```bash
python3 pipeline/scripts/model2/prepare_lamar_handoff.py \
  --input /path/CU5.17_EGFP_GC.treatment_specific.tsv.gz \
  --output /path/CU5.17_EGFP_GC.Lamar_handoff_metadata.tsv
```

With a 101-nt transcript-oriented sequence window:

```bash
python3 pipeline/scripts/model2/prepare_lamar_handoff.py \
  --input /path/CU5.17_EGFP_GC.treatment_specific.tsv.gz \
  --output /path/CU5.17_EGFP_GC.Lamar_handoff_101nt.tsv \
  --reference /path/GRCh38.primary_assembly.genome.fa \
  --flank 50
```

The script verifies that the center of every oriented sequence is `C`.

## Important boundary

Control edit-rate values are missing in the frozen final table because control
sites were not called by REDItools2. Missing calls must not be converted to
zero. Therefore:

- the current table is suitable for **Lamar inference and candidate ranking**;
- `median_treated_edit_rate` is only a provisional training target;
- a background-corrected training label requires control-site base counts and
  depth at the same coordinates;
- fine-tuning also requires sufficiently covered low-editing/background sites;
- train/validation/test splits should be grouped by gene, transcript or genomic
  region to reduce sequence leakage.

The 16 catalogue-overlapping alleles are exclusions, not negative editing
examples.
