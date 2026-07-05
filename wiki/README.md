# Model 1 — RNA-editing evidence pipeline

## Overview

REWIRE recruits a cytidine deaminase to a selected RNA target. Reporter editing establishes on-target activity, but not transcriptome-wide specificity. Model 1 asks:

> **Which C-to-U signals are reproducible after treatment, adequately measured in controls, consistent with transcript orientation and not readily explained by known 293T genomic variation?**

Model 1 produces an auditable set of treatment-associated RNA-editing candidates. It does not claim that every retained mismatch is a biological off-target.

## Design decision

We first evaluated three public WGS BioSamples as a possible genomic blacklist. Only 19.2–26.3% of their reads mapped to GRCh38, strict per-sample filtering retained 137–230 variants, and the two-of-three consensus contained only 118 variants. This was insufficient for genome-wide filtering.

The final workflow therefore uses the `293T_CG` VCF from the [HEK293 Genome Project](https://hek293genome.org/v2/data.php). The database call set is converted from NCBI build 36/hg18 to GRCh38, validated against the same reference used by the RNA branch and joined to RNA candidates by exact allele.

## Input data

| Condition | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

The genomic catalogue is `293T_CG.vcf.gz`, generated with the Complete Genomics pipeline and released in build36/hg18 coordinates.<sup>1</sup>

## Assumptions

1. A mismatch is not automatically an editing event.
2. A missing control call is not negative evidence unless the coordinate has sufficient control depth.
3. Transcript-level C-to-U appears as genomic C→T on positive-strand transcripts and G→A on negative-strand transcripts.
4. The 293T catalogue is external evidence, not WGS matched to the exact CU5.17 cell batch.
5. hg18 and GRCh38 coordinates cannot be compared without assembly harmonization.

## Workflow

![Figure 1. RNA-editing evidence generation pipeline](assets/figure1_model1_evidence_pipeline.svg)

**Figure 1 | RNA-editing evidence generation pipeline.** **a,** RNA-seq produces replicate-, control-, depth- and strand-aware evidence. **b,** The HEK293 Genome Project catalogue is converted from hg18 to GRCh38 and reference-validated. **c,** The branches are joined by `CHROM:POS:REF:ALT` to generate candidates for orthogonal validation.

# Method

## 1. RNA-seq alignment and preprocessing

Paired-end reads are aligned to the GRCh38 primary assembly with STAR in two-pass mode.<sup>2</sup> GATK `MarkDuplicates` records PCR and optical duplicates, and `SplitNCigarReads` processes reads spanning splice junctions.<sup>3</sup> Each final BAM must be sorted, indexed and readable.

## 2. Substitution calling

REDItools2 scans each library independently.<sup>4</sup> The core settings are:

| Setting | Interpretation |
|---|---|
| `-S` | report positions containing a substitution |
| `-me 20` | require at least 20 edited reads at a reported position |
| mapping quality | discard reads below 20 |
| base quality | discard bases below 30 |

The edited-read threshold favours strongly supported events and may miss low-frequency editing.

## 3. Transcript orientation

Union candidates are annotated with Ensembl VEP.<sup>5</sup>

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

Sites assigned to transcripts on both orientations are marked as ambiguous.

## 4. Control evidence

Every union candidate is queried directly in all six BAM files using base quality ≥30 and mapping quality ≥20. This separates **not called despite sufficient depth** from **not observed because the sample was uninformative**.

## 5. 293T catalogue harmonization

**Motivation.** The source catalogue uses hg18 while the RNA branch uses GRCh38.

**Mechanism.** The workflow:

1. retains PASS biallelic SNPs with a non-reference genotype;
2. converts hg18 coordinates to GRCh38 using the UCSC `hg18ToHg38` chain and CrossMap;<sup>6,7</sup>
3. removes unmapped records;
4. validates REF alleles against the project GRCh38 FASTA;
5. removes REF mismatches, normalizes, sorts and tabix-indexes the VCF;
6. extracts a C→T/G→A subset for C-to-U analysis.

**Role.** The catalogue marks exact RNA alleles that are plausible genomic variants. These records remain in the complete site matrix but are excluded from the high-confidence core set.

## 6. Evidence integration

The default high-confidence definition requires:

```text
called in all three treated replicates
AND called in no control replicate
AND candidate-site depth ≥20 in all six RNA-seq libraries
AND transcript orientation consistent with C-to-U editing
AND exact CHROM:POS:REF:ALT absent from the 293T catalogue
```

The final site matrix retains sample-level call status, independent depth, REDItools2 coverage, alternate-read count, editing fraction and catalogue overlap.

## Catalogue quality control

| Processing stage | Variant count |
|---|---:|
| hg18 PASS biallelic SNPs | 2,914,465 |
| CrossMap-unmapped records | 5,979 |
| GRCh38 REF mismatches removed | 22,761 |
| Final GRCh38 PASS biallelic SNPs | 2,885,725 |
| C→T or G→A catalogue alleles | 997,698 |

The final VCF is coordinate-sorted, bgzip-compressed, tabix-indexed and validated against the same GRCh38 FASTA used for the RNA branch.

## Validation and reporting

The workflow checks FASTQ integrity, BAM read groups and indexing, duplicate metrics, REDItools2 interval completeness, depth in all six samples, VEP orientation, liftover success, GRCh38 REF compatibility and exact-allele catalogue overlap.

The final report should provide:

- strand-consistent candidate count;
- treated and control replicate support;
- minimum depth across all six libraries;
- catalogue-overlap count;
- the final high-confidence set;
- the complete site matrix, including excluded records.

## Wet-lab integration

Model 1 prioritizes candidates for targeted amplicon sequencing, independent RNA-seq or another orthogonal assay. Candidates should span multiple editing fractions and sequence contexts rather than only the highest-ranked sites.

## Connection to Model 2 — Lamar

Model 1 provides genomic coordinates, transcript-oriented sequence context, treated and control editing fractions, replicate support, depth and catalogue-overlap status. A continuous label can be defined as:

```python
corrected_efficiency = max(
    0.0,
    median_treated_editing_rate - median_control_editing_rate
)
```

Training and evaluation should be separated by gene, transcript or genomic region rather than by random rows to reduce sequence leakage.

## Limitations

RNA-seq mismatches may reflect genomic variants, alignment ambiguity, sequencing artefacts, endogenous modification, repetitive sequence or batch effects. The `293T_CG` resource represents one database 293T genome and one calling pipeline; it is not matched WGS from the experimental batch. Absence from the catalogue does not prove that a site is free of genomic variation.

Retained sites should be described as:

> **treatment-associated RNA-editing candidates filtered against an external 293T genomic-variant catalogue**

## Contribution

Model 1 combines a fixed three-treated/three-control design, coverage-qualified REDItools2 calling, transcript-oriented interpretation, independent control-depth assessment and reproducible harmonization of a 293T genomic catalogue. The result is an auditable evidence layer for wet-lab prioritization and downstream sequence modelling.

## References

1. Lin, Y.-C. *et al.* Genome dynamics of the human embryonic kidney 293 lineage in response to cell biology manipulations. *Nat. Commun.* **5**, 4767 (2014). doi:10.1038/ncomms5767
2. Dobin, A. *et al.* STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* **29**, 15–21 (2013). doi:10.1093/bioinformatics/bts635
3. McKenna, A. *et al.* The Genome Analysis Toolkit. *Genome Res.* **20**, 1297–1303 (2010). doi:10.1101/gr.107524.110
4. Picardi, E. & Pesole, G. REDItools. *Bioinformatics* **29**, 1813–1814 (2013). doi:10.1093/bioinformatics/btt287
5. McLaren, W. *et al.* Ensembl Variant Effect Predictor. *Genome Biol.* **17**, 122 (2016). doi:10.1186/s13059-016-0974-4
6. Kent, W. J. *et al.* The human genome browser at UCSC. *Genome Res.* **12**, 996–1006 (2002). doi:10.1101/gr.229102
7. Zhao, H. *et al.* CrossMap. *Bioinformatics* **30**, 1006–1007 (2014). doi:10.1093/bioinformatics/btt730
8. Danecek, P. *et al.* SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021). doi:10.1093/gigascience/giab008

**Code and reproducibility:** https://github.com/pdx12320/REWIRE-RNA-editing-pipeline
