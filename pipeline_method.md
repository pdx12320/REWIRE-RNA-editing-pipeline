# RNA-seq-based C-to-U Editing Site Extraction Pipeline

## 1. Why we built this pipeline

The PUF-APOBEC system requires reliable editing labels for downstream prediction models. Raw RNA-seq editing-like signals can contain background sequencing errors, endogenous RNA editing, SNP-like signals, and mapping artifacts. Therefore, we built a treated-control corrected pipeline to identify transcriptome-wide C-to-U candidate editing sites from RNA-seq data.

The pipeline starts from SRA files and produces AI-ready labels for LAMAR-style model training.

## 2. Data processing workflow

```text
SRA accessions
↓
prefetch / fasterq-dump
↓
FASTQ files
↓
STAR alignment to GRCh38
↓
Sorted BAM files
↓
GENCODE exon C-equivalent site scanning
↓
Base counting at each covered site
↓
Replicate count merging
↓
Treated-control correction
↓
Candidate editing sites
↓
Positive and negative sets for LAMAR
```

## 3. Strand-aware C-to-U signal definition

RNA C-to-U editing is represented differently depending on transcript strand because reads are aligned to the genomic reference.

For positive-strand transcripts:

```text
transcript C = genomic C
C-to-U editing appears as genomic C→T
```

For negative-strand transcripts:

```text
transcript C = genomic G
C-to-U editing appears as genomic G→A
```

Thus, the pipeline scans transcript-oriented C-equivalent sites:

```text
+ strand: reference base C, edited base T
- strand: reference base G, edited base A
```

The negative-strand G→A signal is only a genomic-coordinate representation of transcript-level C-to-U editing.

## 4. Site-level editing rate calculation

For each covered C-equivalent site, we calculate:

```text
total_count  = number of reads covering the site
edited_count = number of reads supporting the edited base
edit_rate    = edited_count / total_count
```

Only sites with enough coverage are retained:

```text
total_count >= 20
```

## 5. Replicate merging

Replicates are merged at the count level rather than by averaging editing rates.

```text
merged_edited_count = edited_count_rep1 + edited_count_rep2 + ...
merged_total_count  = total_count_rep1  + total_count_rep2  + ...
merged_edit_rate    = merged_edited_count / merged_total_count
```

This avoids bias caused by unequal coverage across replicates.

## 6. Treated-control correction

The treated and control tables are joined using the same genomic C-equivalent site:

```text
chrom + strand + genomic_pos_1based
```

For each matched site, we calculate:

```text
delta_edit_rate     = treated_edit_rate - control_edit_rate
corrected_edit_rate = (treated_edit_rate - control_edit_rate) / (1 - control_edit_rate)
```

The corrected rate is clipped to the interval [0, 1].

## 7. Candidate editing site definition

A site is called a control-corrected candidate editing site if it passes:

```text
treated_total_count >= 20
control_total_count >= 20
treated_edit_rate >= 0.05
control_edit_rate <= 0.02
delta_edit_rate >= 0.05
corrected_edit_rate >= 0.05
```

For high-confidence candidate sites, we use stricter thresholds:

```text
treated_edit_rate >= 0.10
corrected_edit_rate >= 0.10
control_edit_rate <= 0.02
```

If statistical testing is performed, we further require:

```text
Fisher q <= 0.05
```

## 8. PUF target similarity annotation

Nearby PUF target similarity is recorded as an auxiliary feature, not as a strict filter. For each candidate site, the surrounding sequence can be scanned to report:

```text
best_puf_hamming
best_puf_kmer
best_puf_rel_start
```

Because endogenous off-target candidate sites may not contain a perfect PUF target, PUF mismatch should not be used as a hard cutoff in the first-pass dataset.

## 9. Output files

| File | Description | Downstream use |
|---|---|---|
| `treated.all_C_sites.genomic.cov20.tsv.gz` | All covered C-equivalent sites in treated samples | positive/background site pool |
| `control.all_C_sites.genomic.cov20.tsv.gz` | All covered C-equivalent sites in control samples | background correction |
| `*.edited_sites.genomic.cov20.rate005.tsv.gz` | sites with edit_rate >= 0.05 before control correction | quick exploration |
| `*.strong_sites.genomic.cov20.rate05.tsv.gz` | sites with edit_rate >= 0.5 before control correction | strong case-study sites |
| `*.control_corrected.PUF_candidates.tsv.gz` | treated-control corrected candidates | candidate positive set |

## 10. Model training interpretation

This dataset is better suited for a candidate editing classifier or ranking model than for precise editing-efficiency regression.

Recommended labels:

```text
positive = high-confidence control-corrected candidate sites
negative = clean C sites with treated_edit_rate < 0.02 and control_edit_rate < 0.02
```

The model output should be interpreted as candidate editing risk or ranking score, not as fully validated biochemical editing efficiency.
