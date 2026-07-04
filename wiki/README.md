# Model 1 — RNA-editing evidence pipeline

## Overview

REWIRE recruits a cytidine deaminase to a selected RNA sequence. Reporter editing demonstrates activity at the intended target, but does not establish transcriptome-wide specificity. Model 1 therefore asks a narrower and testable question:

> **Which C-to-U signals are reproducible in treated samples, adequately measured in controls, consistent with transcript orientation and not readily explained by genomic variation?**

Model 1 does not label every RNA-seq mismatch as an off-target. It records an evidence chain for each candidate: read support, replicate consistency, control coverage, transcript orientation and overlap with a public HEK293T genomic-variant catalogue.

---

## Experimental design

We analysed three editor-treated RNA-seq libraries and three control libraries.

| Condition | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

Replicates test whether a signal is reproducible rather than library-specific. Controls provide background evidence only when the same coordinate has sufficient coverage.

Three public HEK293T WGS runs were configured as an external genomic-variant resource:

```text
SRR37832939
SRR37832940
SRR37832941
```

Their metadata are checked before analysis. Runs from one BioSample are merged; runs from different BioSamples are called independently and combined as an exact-allele two-of-three consensus blacklist.

---

## Pipeline overview

![Figure 1. RNA-editing evidence generation pipeline](assets/figure1_model1_evidence_pipeline.svg)

**Figure 1 | RNA-editing evidence generation pipeline.** **a,** The RNA-seq branch identifies quality-supported substitutions and adds transcript orientation, all-sample depth and treated/control evidence. **b,** The WGS branch aligns public HEK293T genomes and constructs a genomic-SNV blacklist. **c,** The branches are integrated by exact `CHROM:POS:REF:ALT` matching to produce treatment-associated C-to-U candidates. Public WGS is used as an external blacklist rather than as WGS matched to the experimental cell batch.

---

## Design principles and assumptions

Model 1 was built around four decision rules.

1. **A mismatch is not automatically an editing event.** RNA-seq mismatches can arise from sequencing error, alignment ambiguity, endogenous editing or genomic variation.
2. **A missing control call is not automatically a negative observation.** The coordinate must be independently covered in the control BAM.
3. **C-to-U editing must be interpreted in transcript orientation.** It appears as genomic C→T on positive-strand transcripts and genomic G→A on negative-strand transcripts.
4. **Public WGS can flag genomic variants but cannot fully replace matched WGS.** HEK293T sublines accumulate different variants during passage and laboratory propagation.

---

# Method

## 1. RNA-seq download and sample tracking

SRA Toolkit downloads each accession and converts it to paired FASTQ files. A fixed manifest stores the sample name, condition, replicate and accession, preventing treated/control labels from being re-entered at later stages.

```bash
python3 pipeline/scripts/rna/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest pipeline/config/samples.tsv \
  --threads 16
```

**Quality checkpoint:** both FASTQ mates must exist and be non-empty.

## 2. Splice-aware alignment with STAR

STAR aligns paired RNA-seq reads to the GRCh38 primary assembly in two-pass mode. The first pass discovers splice junctions; the second pass uses those junctions during final alignment. STAR writes coordinate-sorted BAM files and read-group fields.

```bash
python3 pipeline/scripts/rna/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest pipeline/config/samples.tsv \
  --threads 50
```

Read groups preserve sample, library and platform identity for GATK and Picard. They do not normalise sequencing depth.

## 3. RNA-aware BAM preprocessing

GATK `MarkDuplicates` flags PCR and optical duplicates and records duplicate metrics. `SplitNCigarReads` then processes reads spanning splice junctions into exon-aligned segments suitable for position-level mismatch analysis.

```bash
python3 pipeline/scripts/rna/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest pipeline/config/samples.tsv \
  --java-options=-Xmx16g
```

**Quality checkpoint:** each final BAM must pass `samtools quickcheck` and have a coordinate index.

## 4. Coverage-balanced REDItools2 calling

Parallel REDItools2 uses a per-position coverage map to divide the reference into intervals with similar computational load. This coverage map is a scheduling input, not an editing result. Coverage generation is limited to eight concurrent `samtools depth` jobs, whereas REDItools2 uses 30 MPI processes.

```bash
conda activate reditools2_py2

nohup bash pipeline/scripts/rna/run_reditools_all_samples.sh \
  "$PROJECT" "$REF" "$REDITOOLS" "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

The core REDItools2 settings are:

| Parameter | Interpretation |
|---|---|
| `-S` | report positions containing an observed substitution |
| `-me 20` | require at least 20 edited reads at a reported position |
| default `-q 20` | minimum read mapping quality |
| default `-bq 30` | minimum base quality |
| `-G`, `-D` | complete and per-contig coverage inputs |
| `-np 30` | 30 MPI processes |

`-me 20` is an edited-read threshold, not a total-depth threshold. It favours strongly supported events and may miss low-frequency editing.

### Preserving GRCh38 contig identifiers

GRCh38 contains supplementary contigs such as `GL000194.1`, `GL000205.2` and `KI270750.1`. The version suffix is part of the reference identifier. The original REDItools2 temporary-file parser removed everything after the first dot. We replaced it with:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

This removes only the terminal `.gz` extension and preserves the full contig name.

## 5. Transcript-oriented C-to-U interpretation

REDItools2 reports substitutions in genomic coordinates, whereas editing occurs in RNA transcripts. Candidate substitutions from all six samples are combined into one union VCF and annotated with VEP.

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

The negative-strand G→A representation is the reverse complement of transcript-level C-to-U editing; it does not imply biochemical editing of G. Coordinates assigned to transcripts on both orientations are treated as ambiguous.

## 6. Independent depth and control subtraction

A coordinate absent from a control REDItools2 table may simply have insufficient coverage. We therefore query every union candidate in every treated and control BAM using base quality ≥30 and mapping quality ≥20.

```bash
bash pipeline/scripts/rna/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  pipeline/config/samples.tsv
```

This distinguishes **not called despite sufficient depth** from **not observed because the sample was uninformative**.

## 7. Public HEK293T WGS filtering

The WGS branch uses SRA Toolkit, BWA-MEM2 or BWA-MEM, GATK `MarkDuplicates` and bcftools `mpileup/call`. All variants are normalised against the same GRCh38 FASTA used for RNA-seq.

```bash
conda env create -f pipeline/env/wgs_pipeline.yml
conda activate rewire_wgs

nohup bash pipeline/scripts/wgs/run_3run_wgs_pipeline.sh \
  --runs pipeline/config/wgs_runs.tsv \
  --reference "$REF" \
  --outdir "$WGS_OUT" \
  --mode auto \
  --threads 32 \
  --min-dp 10 \
  --min-alt 3 \
  --min-vaf 0.05 \
  --min-qual 20 \
  > "$WGS_OUT.pipeline.log" 2>&1 &
```

A single-run WGS SNV is retained when depth is ≥10, alternate-read count is ≥3, alternate-allele fraction is ≥0.05 and QUAL is ≥20. If runs represent different BioSamples, the conservative blacklist retains exact alleles found in at least two of three call sets.

## 8. Evidence integration

The final comparison joins RNA and WGS evidence using the exact chromosome, position, reference and alternate allele.

```bash
python3 pipeline/scripts/rna/filter_c_to_u_and_compare.py \
  --manifest pipeline/config/samples.tsv \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --wgs-vcf "$WGS_BLACKLIST" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

The conservative default definition requires:

```text
called in all three treated replicates
AND called in no control replicate
AND candidate-site depth ≥20 in all six samples
AND transcript orientation consistent with C-to-U editing
AND the exact allele absent from the selected WGS blacklist
```

The evidence matrix retains call status, independent depth, REDItools2 depth, alternate-read count and editing fraction for every sample. Thresholds can therefore be re-evaluated without repeating alignment or REDItools2 discovery.

---

## Validation and quality control

The workflow includes explicit checkpoints at each boundary:

- FASTQ mate integrity after SRA conversion;
- read-group presence and coordinate sorting after STAR;
- duplicate metrics and BAM integrity after GATK;
- complete coverage and interval files before REDItools2 merging;
- preservation of versioned GRCh38 contig names;
- valid bgzip and tabix outputs;
- informative depth at each candidate in all six RNA-seq samples;
- exact-allele rather than position-only WGS overlap.

These checks do not prove biological editing, but they make technical failures visible and the candidate-selection logic auditable.

---

## Results

Numerical results are not included in this repository until all six RNA-seq samples and the selected WGS workflow complete the same quality-control procedure. The final wiki will report:

1. REDItools2 calls per sample;
2. strand-consistent C-to-U candidates;
3. treated-replicate overlap;
4. control depth and alternate-read support;
5. candidates removed by the public HEK293T blacklist;
6. the final ranked candidate set.

---

## Integration with the wet lab

Model 1 does not replace experimental validation. It reduces the candidate space and records why each site was prioritised. High-ranking candidates can be tested by targeted amplicon sequencing, independent RNA-seq or another orthogonal assay. The same evidence table can also provide carefully defined positive and background examples for downstream sequence models.

---

## Limitations

RNA-seq mismatches may arise from genomic variants, alignment ambiguity, sequencing artefacts, endogenous RNA modification, repetitive sequence or batch effects. Treated/control comparison, independent depth and public WGS filtering reduce these alternatives but do not eliminate them.

Because the WGS data are public and were not generated from the exact CU5.17 experimental cell batch, retained sites should be described as:

> **treatment-associated RNA-editing candidates filtered against an external HEK293T genomic-variant catalogue**

They should not be described as definitively SNV-free off-targets. High-priority sites still require orthogonal validation.

---

## Contribution

Model 1 provides:

- a fixed three-treated/three-control analysis design;
- splice-aware alignment and RNA-specific BAM preprocessing;
- coverage-balanced MPI REDItools2 calling;
- a fix for versioned GRCh38 contig parsing;
- transcript-oriented C-to-U interpretation;
- independent depth assessment in all six samples;
- treated/control evidence integration;
- an automated three-run public HEK293T WGS workflow;
- an exact-allele merged or two-of-three SNV blacklist;
- an auditable site-level evidence matrix.

---

## References

1. Dobin A. *et al.* STAR: ultrafast universal RNA-seq aligner. **Bioinformatics** 29, 15–21 (2013).
2. McKenna A. *et al.* The Genome Analysis Toolkit: a MapReduce framework for analysing next-generation DNA sequencing data. **Genome Research** 20, 1297–1303 (2010).
3. Picardi E. & Pesole G. REDItools: high-throughput RNA editing detection made easy. **Bioinformatics** 29, 1813–1814 (2013).
4. McLaren W. *et al.* The Ensembl Variant Effect Predictor. **Genome Biology** 17, 122 (2016).
5. Li H. & Durbin R. Fast and accurate short read alignment with Burrows–Wheeler transform. **Bioinformatics** 25, 1754–1760 (2009).
6. Danecek P. *et al.* Twelve years of SAMtools and BCFtools. **GigaScience** 10, giab008 (2021).

---

**Code and reproducibility.** Complete scripts, manifests, environment files, expected outputs and troubleshooting notes are available at:  
**https://github.com/pdx12320/REWIRE-RNA-editing-pipeline**
