# Dry Lab — Model 1: RNA-editing evidence pipeline

## Motivation

The REWIRE editor is designed to recruit a cytidine deaminase to a selected RNA sequence. Reporter editing demonstrates activity at the intended site, but transcriptome-wide specificity requires a separate analysis. Model 1 therefore converts matched treated and control RNA-seq data into a replicate-aware matrix of candidate C-to-U events.

## Evidence chain

```text
raw paired-end RNA-seq
→ alignment and RNA-specific preprocessing
→ quality-filtered substitution discovery
→ transcript-strand annotation
→ independent depth confirmation
→ treated/control replicate comparison
→ optional genomic-variant removal
```

The design deliberately separates a raw RNA-seq mismatch from a final candidate. A retained site must be supported by edited reads, interpreted correctly in transcript orientation, reproduced across treated libraries, and evaluated in controls at the same coordinate.

## Why this is Model 1

Model 1 generates the evidence used by the rest of the dry-lab framework. It provides traceable candidate labels rather than relying on an unfiltered mismatch table. This is important for downstream prediction because model performance is only meaningful when the training labels have a defensible biological origin.

The full Wiki-ready page is available in `Model1_RNA_editing_evidence_pipeline.md`, and all reproducibility commands are listed in `Code_and_commands.md`.
