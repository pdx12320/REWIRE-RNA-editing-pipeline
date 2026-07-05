# Model 1 — RNA-editing evidence pipeline

## Overview

REWIRE recruits a cytidine deaminase to a selected RNA target. Reporter editing establishes on-target activity, but not transcriptome-wide specificity. Model 1 asks:

> **Which C-to-U signals are reproducible after treatment, absent from control calls, consistent with transcript orientation and not readily explained by known 293T genomic variation?**

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

## Assumptions and boundaries

1. A mismatch is not automatically an editing event.
2. A missing control call is not equivalent to zero editing or adequate control coverage.
3. Transcript-level C-to-U appears as genomic C→T on positive-strand transcripts and G→A on negative-strand transcripts.
4. The 293T catalogue is external evidence, not WGS matched to the exact CU5.17 cell batch.
5. hg18 and GRCh38 coordinates cannot be compared without assembly harmonization.

## Workflow

![Figure 1. RNA-editing evidence generation pipeline](assets/figure1_model1_evidence_pipeline.svg)

**Figure 1 | RNA-editing evidence generation pipeline.** **a,** RNA-seq produces replicate-, control-call- and strand-aware evidence; all-sample candidate depth is a recommended validation layer. **b,** The HEK293 Genome Project catalogue is converted from hg18 to GRCh38 and reference-validated. **c,** The branches are joined by `CHROM:POS:REF:ALT` to generate a catalogue-filtered screening set for downstream ranking and orthogonal validation.

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

## 4. Treated/control call comparison

The frozen screening table retains sites called in all three treated replicates and not called in any of the three controls under the original REDItools2 filtering settings. Because the final legacy output does not contain independent candidate-site depth and base counts for the uncalled controls, these records are treated as **screening candidates**, not as proof of control absence.

A stricter future gate should query each candidate directly in all six BAM files. This distinguishes **not called despite sufficient depth** from **not observed because the sample was uninformative**.

## 5. 293T catalogue harmonization

**Motivation.** The source catalogue uses hg18 while the RNA branch uses GRCh38.

**Mechanism.** The workflow:

1. retains PASS biallelic SNPs with a non-reference genotype;
2. converts hg18 coordinates to GRCh38 using the UCSC `hg18ToHg38` chain and CrossMap;<sup>6,7</sup>
3. removes unmapped records;
4. validates REF alleles against the project GRCh38 FASTA;
5. removes REF mismatches, normalizes, sorts and tabix-indexes the VCF;
6. performs exact `CHROM:POS:REF:ALT` comparison with RNA candidates.

**Role.** The catalogue marks RNA alleles that are plausible genomic variants. These records are retained in an exclusion table rather than silently discarded.

## 6. Exact-allele integration

For the frozen result set, the implemented screening definition is:

```text
called in all three treated replicates
AND called in no control replicate under the original REDItools2 filter
AND transcript orientation consistent with C-to-U editing
AND exact CHROM:POS:REF:ALT absent from the 293T catalogue
```

This definition should not be conflated with a fully depth-qualified high-confidence off-target set.

## Catalogue quality control

| Processing stage | Variant count |
|---|---:|
| hg18 PASS biallelic SNPs | 2,914,465 |
| CrossMap-unmapped records | 5,979 |
| GRCh38 REF mismatches removed | 22,761 |
| Final GRCh38 PASS biallelic SNPs | 2,885,725 |

The final VCF is coordinate-sorted, bgzip-compressed, tabix-indexed and validated against the same GRCh38 FASTA used for the RNA branch.

# Results

The evidence funnel produced:

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact alleles overlapping the GRCh38-harmonized 293T catalogue | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

Thus, exact catalogue comparison removed 16 of the 3,349 treatment-specific candidates, leaving 3,333 candidates for sequence-context analysis, Lamar inference and experimental prioritization.

The final output is not described as a definitive off-target list because independent depth and base counts for the uncalled controls are absent from the frozen table.

## Validation and reporting

The implemented workflow checks FASTQ integrity, BAM read groups and indexing, duplicate metrics, REDItools2 interval completeness, VEP orientation, liftover success, GRCh38 REF compatibility and exact-allele catalogue overlap.

The complete report should preserve:

- the 3,333 retained screening candidates;
- the 16 catalogue-overlapping exclusions;
- treated replicate coverage, alternate-read count and editing rate;
- strand orientation and exact genomic allele;
- the limitation that control non-calls lack independent depth evidence.

## Wet-lab integration

Model 1 prioritizes candidates for targeted amplicon sequencing, independent RNA-seq or another orthogonal assay. Candidates should span multiple editing fractions and sequence contexts rather than only the highest-ranked sites.

## Connection to Model 2 — Lamar

The 3,333 retained records can be used immediately for Lamar **inference and ranking** after fixed-length transcript-oriented sequence extraction. Positive-strand C→T sites retain the genomic sequence orientation; negative-strand G→A sites must be reverse-complemented so that the center is transcript-level C.

For this frozen table, `median_treated_edit_rate` is the recommended provisional continuous target because it is robust across the three treated replicates. However, control edit rates are missing because the control sites were not called. Missing values must not be converted to zero.

Therefore, model fine-tuning requires additional data:

- direct control-site base counts and depth at the same coordinates;
- sufficiently covered low-editing or unedited background sites;
- gene-, transcript- or region-grouped train/validation/test splits to reduce sequence leakage.

The 16 catalogue-overlapping alleles are exclusions, not negative editing examples.

## Limitations

RNA-seq mismatches may reflect genomic variants, alignment ambiguity, sequencing artefacts, endogenous modification, repetitive sequence or batch effects. The `293T_CG` resource represents one database 293T genome and one calling pipeline; it is not matched WGS from the experimental batch. Absence from the catalogue does not prove that a site is free of genomic variation.

Retained sites should be described as:

> **catalogue-filtered treatment-associated C-to-U screening candidates**

They should not be described as definitively SNV-free off-targets or fully depth-qualified editing events.

## Contribution

Model 1 combines a fixed three-treated/three-control call design, transcript-oriented interpretation, reproducible harmonization of a 293T genomic catalogue and exact-allele integration. The result is an auditable screening layer for wet-lab prioritization and downstream Lamar inference, with the remaining control-depth uncertainty stated explicitly.

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
