# Transcriptome-wide C-to-U RNA-editing evidence pipeline

## Overview

REWIRE directs a PUF–APOBEC editor to a selected RNA target. Reporter editing confirms on-target activity, but it does not reveal transcriptome-wide off-target signals. We therefore built a screening pipeline to identify C-to-U candidates that are reproducible after treatment, absent from control calls, consistent with transcript strand and not already present in a 293T genomic-variant catalogue.

The final output is an auditable screening set. It is not a definitive off-target list.

## Final workflow

![RNA-editing evidence generation pipeline](assets/figure1_model1_evidence_pipeline.svg)

1. Align three treated and three control RNA-seq libraries to GRCh38.
2. Preprocess the BAM files and call substitutions with REDItools2.
3. Use VEP to interpret genomic C→T and G→A events in transcript orientation.
4. Retain sites called in all three treated replicates and not called in the controls under the original filter.
5. Remove exact `CHROM:POS:REF:ALT` matches to the GRCh38-harmonized `293T_CG` catalogue.

The full commands, failure logs and implementation decisions are documented in the [engineering notebook](../docs/ENGINEERING_CYCLE.md) and [pipeline guide](../pipeline/README.md).

## Engineering cycle

### Cycle 1 — Establish the RNA evidence branch

**Design.** Analyse each replicate independently so that reproducibility remains visible.

**Build.** We combined STAR, GATK, REDItools2 and VEP to generate strand-consistent candidate alleles.

**Test.** The pipeline produced 9,930 strand-consistent sites. Of these, 4,778 were called in all three treated replicates, and 3,349 were not called in any control under the original REDItools2 filter.

**Learn.** A control non-call is not the same as zero editing. The retained sites therefore remain screening candidates unless control-site depth and base counts are measured directly.

### Cycle 2 — Test public WGS as a blacklist

**Design.** We first attempted to reconstruct a 293T genomic blacklist from three public WGS runs.

**Build.** The runs were aligned and filtered separately, then compared by exact allele.

**Test.** Only 19.2–26.3% of reads mapped to GRCh38. The two-of-three consensus contained 118 variants.

**Learn.** This catalogue was too sparse for genome-wide exclusion, so the public-WGS route was removed from the final workflow.

### Cycle 3 — Harmonize the 293T_CG catalogue

**Design.** Replace the failed WGS reconstruction with the database-released `293T_CG` VCF from the HEK293 Genome Project.<sup>1</sup>

**Build.** The source VCF was filtered, lifted from hg18 to GRCh38 with CrossMap, checked against the project FASTA, normalized, sorted and indexed.

**Test.** The conversion retained 2,885,725 GRCh38 SNPs. During development, we also detected a compressed VCF with a misleading extension, a failed chain download, unsorted CrossMap output and 22,761 REF mismatches.

**Learn.** The final script detects compression from file content, accepts a local chain file, validates REF alleles and sorts before indexing.

### Cycle 4 — Integrate and audit

**Design.** Compare RNA candidates with the catalogue by exact allele and keep excluded records for audit.

**Build.** A compatibility filter was added for the completed legacy treatment-specific table, which lacked the newer depth-summary column.

**Test.** Sixteen of 3,349 treatment-specific candidates matched the 293T catalogue, leaving 3,333 retained sites.

**Learn.** Catalogue overlap supports exclusion as a plausible genomic variant, but catalogue absence does not prove that a site is SNV-free in the exact experimental subline.

## Results

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact 293T catalogue overlaps | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

The 16 catalogue-overlapping alleles were retained in a separate exclusion table. The 3,333 remaining sites form the final screening set for sequence-context analysis and experimental prioritization.

## Contribution

This work provides:

- a reproducible RNA-editing evidence pipeline for three treated and three control RNA-seq libraries;
- strand-aware interpretation of C-to-U events;
- a tested hg18-to-GRCh38 conversion workflow for the `293T_CG` catalogue;
- exact-allele integration with explicit retention of excluded records;
- a documented DBTL history showing which strategies failed, why they failed and how the workflow changed.

All scripts, commands, QC counts and troubleshooting notes are available in this repository.

## Limitations

The final legacy table does not contain independent depth and base counts for control non-calls. The 3,333 retained records are therefore screening candidates rather than fully depth-qualified editing events.

The `293T_CG` catalogue is an external 293T resource, not matched WGS from the exact experimental batch. A retained site may still represent a subline-specific genomic variant.

RNA-seq mismatches can also arise from alignment ambiguity, sequencing artefacts, repetitive sequence or endogenous RNA modification. Orthogonal validation is still required.

## References

1. Lin, Y.-C. *et al.* Genome dynamics of the human embryonic kidney 293 lineage in response to cell biology manipulations. *Nat. Commun.* **5**, 4767 (2014).
2. Dobin, A. *et al.* STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* **29**, 15–21 (2013).
3. McKenna, A. *et al.* The Genome Analysis Toolkit. *Genome Res.* **20**, 1297–1303 (2010).
4. Picardi, E. & Pesole, G. REDItools. *Bioinformatics* **29**, 1813–1814 (2013).
5. McLaren, W. *et al.* The Ensembl Variant Effect Predictor. *Genome Biol.* **17**, 122 (2016).
6. Zhao, H. *et al.* CrossMap. *Bioinformatics* **30**, 1006–1007 (2014).
7. Danecek, P. *et al.* Twelve years of SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021).
