# LAMAR binary label design

## Purpose

This document records the model-facing label mechanism used by the current LAMAR binary classification release. It supersedes the earlier called-site funnel as the primary dataset description.

The task predicts whether a 101-nucleotide sequence context should enter experimental validation. It does not predict editing efficiency or establish experimental editability.

## Inputs

- three treated and three control MarkDuplicates BAMs;
- GRCh38 primary-assembly sequence;
- GENCODE v50 primary-assembly exon annotation;
- the broad called-site matrix used only as the positive candidate source and ambiguous-site exclusion set;
- two HEK293T WGS variant resources used for central-position exclusion.

All six BAMs were counted with one pileup implementation. The core filters were mapping quality at least 30, base quality at least 20, duplicate exclusion, disallowed-alignment exclusion, and paired-read overlap control.

## Shared sequence and variant checks

Both classes required a complete transcript-oriented 101-nucleotide context with C at zero-based index 50. Sites with ambiguous orientation, N-containing sequence, invalid center base, or duplicate genomic key were excluded.

A site was excluded when its genomic center occurred in either WGS resource, regardless of the reported alternative allele.

Low complexity was defined by an OR rule. A site failed when the 101-nucleotide context met any of these conditions:

- base-2 single-nucleotide Shannon entropy below 1.20;
- longest homopolymer run at least 20 nucleotides;
- maximum deterministic dinucleotide-repeat coverage at least 0.80 across either repeat phase.

External basewise mappability was unavailable and was recorded as `NA_RESOURCE_MISSING`. Low-complexity filtering is not an equivalent substitute for external mappability validation.

## Computational positives

Positive candidates were recounted from the broad called-site matrix. A site entered `positive_main` when it satisfied every condition below:

1. corrected editing efficiency was strictly greater than 0.10;
2. at least two of three treated replicates had usable depth at least 20;
3. at least two of three control replicates had usable depth at least 20;
4. the median control editing rate was at most 0.02;
5. the shared sequence, complexity, orientation, and WGS checks passed.

The frozen release contained 1,513 computational positives. The 1,457-site high-confidence subset additionally required complete six-sample coverage, screening BH-FDR below 0.05, treated MAD at most 0.05, and control MAD at most 0.02.

The screening p-value and BH-FDR are pooled-read screening statistics. They are not replicate-level experimental validation.

## Strict computational negatives

Negative candidates were enumerated from transcript-oriented C centers in the GENCODE v50 exon union. A site entered the strict negative universe only when it satisfied every condition below:

1. usable depth was at least 20 in all six samples;
2. target-ALT count was exactly zero in all six samples;
3. the site did not occur in `positive_main`, the positive sensitivity set, or the broad called-site matrix;
4. the shared sequence, complexity, orientation, and WGS checks passed.

The frozen strict universe contained 2,821,734 sites. These are expressed and observed strict computational negatives, not experimentally proven non-editable positions.

## Leakage-aware splitting

Sites were grouped before splitting. A shared gene, genomic center, overlapping same-strand window, or exact 101-nucleotide sequence forced examples into one leakage group.

The fixed release contained 1,028 train positives, 159 development positives, 165 calibration positives, and 161 held-out positives. The public release audited zero cross-split overlap by gene, genomic key, exact sequence, and leakage group.

## Model-facing boundary

The default model input is `sequence_context` only. Coverage, expression proxies, WGS status, annotation, and other metadata remain available for filtering, matching, stratification, and audit.

The canonical frozen artifacts and complete manifests are maintained in the [PUF fine-tuning repository](https://github.com/pdx12320/PUF_fine_tuning_version1), source commit `77a2a02`.
