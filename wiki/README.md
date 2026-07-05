# Dry Lab — Transcriptome-wide C-to-U RNA-editing screening

## Overview

REWIRE uses a programmable PUF–APOBEC editor to modify a selected RNA target. Reporter assays can confirm activity at the intended site, but they cannot determine whether similar C-to-U changes occur elsewhere in the transcriptome.

We therefore developed a dry-lab pipeline that compares three treated and three control RNA-seq libraries, interprets substitutions in transcript orientation, and removes exact matches to an external 293T genomic-variant catalogue. The output is a ranked, auditable set of treatment-associated C-to-U screening candidates for downstream validation.

## Aim

Our pipeline was designed to answer one question:

> Which transcriptome-wide C-to-U signals are reproducibly detected after treatment, absent from control calls under the same settings, consistent with transcript strand, and not readily explained by known 293T genomic variation?

The analysis does not classify every retained mismatch as a confirmed off-target. Instead, it narrows a large RNA-seq call set into a smaller evidence-based screening set.

## Input data

### RNA-seq libraries

| Condition | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

Each library was processed independently so that replicate support remained visible throughout the analysis.

### Genomic-variant catalogue

We used the `293T_CG` variant catalogue released by the HEK293 Genome Project.<sup>1</sup> The source VCF was generated on NCBI build 36/hg18 coordinates, whereas the RNA-seq branch used GRCh38. We therefore converted and validated the catalogue before integration.

## Workflow

![RNA-editing evidence generation pipeline](assets/figure1_model1_evidence_pipeline.svg)

The pipeline contains two evidence branches:

```text
RNA-seq evidence
    alignment → preprocessing → substitution calling
    → transcript orientation → treated/control comparison

293T genomic evidence
    source filtering → hg18-to-GRCh38 conversion
    → REF validation → exact-allele catalogue

Final integration
    RNA candidate CHROM:POS:REF:ALT
    compared with catalogue CHROM:POS:REF:ALT
```

## 1. RNA-seq alignment and preprocessing

Paired-end reads were aligned to the GRCh38 primary assembly with STAR in two-pass mode.<sup>2</sup> GATK `MarkDuplicates` was used to record duplicate reads, and `SplitNCigarReads` processed reads spanning splice junctions.<sup>3</sup>

The same GRCh38 reference was used throughout the RNA and catalogue branches to avoid assembly-dependent coordinate and allele inconsistencies.

## 2. Substitution calling

REDItools2 was run independently for all six libraries.<sup>4</sup>

Core settings were:

| Parameter | Function |
|---|---|
| `-S` | report positions containing substitutions |
| `-me 20` | require at least 20 edited reads for a reported call |
| mapping quality ≥20 | remove poorly mapped reads |
| base quality ≥30 | remove low-confidence bases |

These settings prioritised strongly supported calls. They may reduce sensitivity to low-frequency editing, but they limit the initial candidate set to sites with substantial read support.

## 3. Transcript-oriented C-to-U interpretation

RNA editing must be interpreted relative to transcript strand rather than genomic substitution alone. VEP was used to annotate candidate alleles and transcript orientation.<sup>5</sup>

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

Candidates inconsistent with transcript-level C-to-U editing were removed. Sites assigned to conflicting transcript orientations were treated as ambiguous rather than forced into one class.

## 4. Replicate and control filtering

A candidate was retained in the treatment-specific screening table when it was:

```text
called in all three treated replicates
AND not called in any of the three controls
AND consistent with transcript-level C-to-U editing
```

The final tables retain replicate-level coverage, alternate-read count and editing rate. This allows later prioritisation to consider both editing magnitude and agreement among treated replicates.

A control non-call is not equivalent to confirmed zero editing. REDItools2 may omit lower-level events that do not reach the edited-read threshold. This point is included in the limitations below.

## 5. 293T catalogue harmonisation

The source `293T_CG` VCF was processed as follows:

```text
retain PASS biallelic SNPs
→ convert hg18 coordinates to GRCh38 with CrossMap
→ remove unmapped records
→ validate REF alleles against the project GRCh38 FASTA
→ remove REF mismatches
→ normalize and coordinate-sort
→ bgzip-compress and tabix-index
```

Catalogue quality control produced:

| Processing stage | Variant count |
|---|---:|
| Source PASS biallelic SNPs | 2,914,465 |
| CrossMap-unmapped records | 5,979 |
| GRCh38 REF mismatches removed | 22,761 |
| Final GRCh38 PASS biallelic SNPs | 2,885,725 |

The final catalogue is an external 293T genomic resource. It is not whole-genome sequencing matched to the exact CU5.17 experimental cell batch.

## 6. Exact-allele integration

RNA candidates were compared with the harmonised catalogue using exact:

```text
CHROM : POS : REF : ALT
```

Coordinate-only matching was not used because different alternate alleles can occur at the same genomic position.

Catalogue-overlapping records were written to a separate exclusion table rather than deleted silently. This preserves a complete record of which candidates were removed and why.

## Results

The successive evidence filters reduced the call set as follows:

| Evidence layer | Number of sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact 293T catalogue overlaps | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

The catalogue comparison removed 16 exact alleles from the 3,349 treatment-specific candidates, leaving 3,333 retained sites.

The final result is provided together with:

- the full annotated pre-catalogue treatment-specific table;
- the 16 catalogue-overlapping exclusions;
- the 3,333 retained screening candidates;
- a summary table recording the filtering counts.

## Interpretation and use

The 3,333 retained sites represent transcriptome-wide candidates that passed the implemented treated-replicate, control-call, transcript-orientation and genomic-catalogue filters.

They can be used to:

- prioritise loci for targeted amplicon sequencing;
- select candidates spanning different editing rates and sequence contexts;
- compare predicted and experimentally measured editing behaviour;
- guide future matched-DNA and independent RNA validation.

They should not be interpreted as 3,333 confirmed biological off-targets.

## Contribution

Our dry-lab work contributes:

1. **A replicate-aware RNA-editing workflow.** Three treated and three control libraries are analysed independently rather than pooled.
2. **Transcript-strand-aware interpretation.** Both genomic C→T and G→A events are correctly mapped to transcript-level C-to-U editing.
3. **A reproducible 293T catalogue conversion.** The public hg18 catalogue is converted to GRCh38, reference-validated, normalised and indexed.
4. **Exact-allele genomic filtering.** RNA candidates are compared using `CHROM:POS:REF:ALT`, reducing coordinate-only false matches.
5. **Auditable outputs.** Retained and excluded sites are reported separately, together with QC counts and reproducible scripts.

The complete implementation, development records and troubleshooting notes are available in the project repository.

## Limitations

### Control evidence

The frozen legacy result table does not contain independent candidate-site depth and base counts for control non-calls. A site absent from the control REDItools2 table may still contain lower-level alternate reads. The retained sites are therefore described as screening candidates rather than fully depth-qualified editing events.

### Genomic evidence

The `293T_CG` catalogue comes from an external 293T genome and calling pipeline. An exact overlap supports exclusion as a plausible genomic variant, but absence from the catalogue does not prove that the exact experimental subline lacks that variant.

### RNA-seq evidence

RNA-seq mismatches may arise from alignment ambiguity, sequencing artefacts, repetitive regions, endogenous RNA modification or batch effects. Targeted amplicon sequencing, matched DNA or an independent RNA-seq experiment is required for orthogonal confirmation.

### Sensitivity

The stringent edited-read threshold favours strong events and may miss low-frequency editing. The final candidate set therefore emphasises specificity over complete sensitivity.

## Reproducibility

The repository separates the public-facing analysis from the full technical record:

- `pipeline/` contains executable scripts and commands;
- `dbtl/` records iterative development, failed approaches and decisions;
- `results/` contains the frozen filtering summary;
- `pipeline/CATALOGUE_PROVENANCE.md` records catalogue source and QC.

## References

1. Lin, Y.-C. *et al.* Genome dynamics of the human embryonic kidney 293 lineage in response to cell biology manipulations. *Nature Communications* **5**, 4767 (2014).
2. Dobin, A. *et al.* STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* **29**, 15–21 (2013).
3. McKenna, A. *et al.* The Genome Analysis Toolkit. *Genome Research* **20**, 1297–1303 (2010).
4. Picardi, E. & Pesole, G. REDItools. *Bioinformatics* **29**, 1813–1814 (2013).
5. McLaren, W. *et al.* The Ensembl Variant Effect Predictor. *Genome Biology* **17**, 122 (2016).
6. Zhao, H. *et al.* CrossMap: a versatile tool for coordinate conversion between genome assemblies. *Bioinformatics* **30**, 1006–1007 (2014).
7. Danecek, P. *et al.* Twelve years of SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021).
