# Interpretation and limitations

## What the pipeline supports

Model 1 identifies RNA-seq positions with evidence consistent with transcript-level C-to-U editing. The strongest candidates are reproducible across treated replicates, absent from matched controls with adequate coverage, strand-consistent after VEP annotation, and absent from the optional matched genomic-variant set.

## What the pipeline does not prove

A computational candidate is not automatically a confirmed biological off-target. RNA-seq mismatches can also arise from:

- genomic variants;
- alignment ambiguity;
- sequencing error;
- endogenous RNA modification;
- repetitive or low-complexity sequence;
- library-specific artifacts;
- sample-specific expression differences.

Matched controls, replicate consistency, depth confirmation, and WGS filtering reduce these alternatives but do not eliminate them completely.

## Sensitivity limitation

The REDItools2 setting `-me 20` requires at least 20 edited reads. This favors well-covered and relatively high-frequency events. Genuine low-frequency editing may therefore be missed. A later sensitivity analysis can use lower edited-read thresholds, but it should also apply stronger artifact controls and independent validation.

## Validation

High-priority candidates should be confirmed with an orthogonal assay. Suitable options include targeted amplicon sequencing, independent RNA-seq, Sanger sequencing for high-frequency events, or another site-specific validation method.

Accordingly, Wiki and presentation text should use the phrase **computational RNA-editing candidates** unless independent experimental validation has been completed.
