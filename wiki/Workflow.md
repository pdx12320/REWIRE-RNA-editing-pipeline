# Computational workflow

## Input

Model 1 begins with three editor-treated paired-end RNA-seq libraries and three matched controls. All samples are processed with the same reference, software, thresholds, and file-integrity checks.

## Processing sequence

```text
SRA
→ FASTQ
→ STAR two-pass alignment
→ GATK MarkDuplicates
→ GATK SplitNCigarReads
→ REDItools2 coverage map
→ parallel substitution calling
→ union VCF and BED
→ VEP transcript-strand annotation
→ all-replicate depth measurement
→ treated/control evidence matrix
```

## Design principle

The workflow separates discovery from interpretation. REDItools2 first reports quality-supported substitutions without assuming that every substitution is caused by the engineered editor. VEP then provides transcript orientation, independent depth checks establish whether each sample was informative at the site, and replicate-level filtering determines whether the signal is treatment-associated.

This staged design avoids three common errors:

1. treating every treated-sample mismatch as editor activity;
2. treating an uncovered control position as a negative observation;
3. ignoring the reverse-complement representation of C-to-U editing on negative-strand transcripts.
