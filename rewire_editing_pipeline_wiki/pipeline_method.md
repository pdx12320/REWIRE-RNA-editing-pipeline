# Pipeline Method

## 1. Input data

The pipeline uses paired RNA-seq samples from an editor-treated group and a matched mock/control group.

Input files:

```text
SRA accessions
Reference genome: GRCh38 primary assembly
Annotation: GENCODE primary assembly GTF
```

## 2. SRA to BAM processing

SRA files are converted to FASTQ files using `prefetch` and `fasterq-dump`. FASTQ files are aligned to GRCh38 using STAR. The output is a coordinate-sorted BAM file for each sample.

```text
SRA → FASTQ → STAR sorted BAM → BAM index
```

For downstream compatibility, STAR BAM files can be linked as `*.split.bam`. In this lightweight version, the `*.split.bam` file is a naming-compatible alias rather than a GATK SplitNCigarReads output.

## 3. Transcript-oriented C-equivalent site scanning

The pipeline scans merged GENCODE exon intervals and identifies transcript-oriented C-equivalent sites.

```text
+ strand: genomic reference base C
- strand: genomic reference base G
```

At each site, aligned reads are counted to calculate the number of reads supporting the reference base and the edited base.

## 4. Editing signal definition

For positive-strand transcripts:

```text
transcript C-to-U editing = genomic C→T signal
```

For negative-strand transcripts:

```text
transcript C-to-U editing = genomic G→A signal
```

The negative-strand G→A signal is not G editing. It is the reverse-complement representation of transcript-level C-to-U editing in genomic coordinates.

## 5. Site-level editing rate calculation

For each covered site:

```text
total_count = number of reads covering the site
edited_count = number of reads carrying the edited base
edit_rate = edited_count / total_count
```

Sites are retained if:

```text
total_count >= 20
```

## 6. Replicate merging

Replicates are merged at the count level:

```text
merged_edited_count = sum(edited_count across replicates)
merged_total_count = sum(total_count across replicates)
merged_edit_rate = merged_edited_count / merged_total_count
```

Editing rates are not averaged directly across replicates, because coverage can differ between samples.

## 7. Treated-control correction

The treated and control tables are joined by genomic coordinate:

```text
chrom + strand + genomic_pos_1based
```

For each matched site:

```text
delta_edit_rate = treated_edit_rate - control_edit_rate
corrected_edit_rate = (treated_edit_rate - control_edit_rate) / (1 - control_edit_rate)
```

## 8. Candidate site filtering

A candidate editing site is defined by:

```text
treated_total_count >= 20
control_total_count >= 20
treated_edit_rate >= 0.05
control_edit_rate <= 0.02
delta_edit_rate >= 0.05
corrected_edit_rate >= 0.05
```

A stricter high-confidence subset can be defined by:

```text
treated_edit_rate >= 0.10
corrected_edit_rate >= 0.10
control_edit_rate <= 0.02
Fisher q <= 0.05, if available
```

## 9. PUF target similarity annotation

The sequence surrounding each candidate site can be scanned for the closest PUF target-like k-mer. The following fields can be recorded:

```text
best_puf_hamming
best_puf_kmer
best_puf_rel_start
```

PUF target similarity is used as an annotation rather than a hard filtering criterion.

## 10. Training data construction

Candidate or high-confidence sites are used as positive samples.

Clean negative sites are sampled from the all-C background table:

```text
treated_edit_rate < 0.02
control_edit_rate < 0.02
total_count >= 20
```

The resulting dataset supports candidate editing classification or genome-wide ranking of potential C-to-U editing sites.
