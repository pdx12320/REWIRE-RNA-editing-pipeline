# Model 1 figure and table plan

This file contains Wiki-ready captions and layouts. The method figures can be prepared before the final numerical results are available.

## Figure 1 — RNA-editing evidence workflow

### Panels

1. Six paired-end RNA-seq libraries: three treated and three controls.
2. STAR alignment to GRCh38.
3. GATK duplicate marking and SplitNCigarReads.
4. REDItools2 coverage-aware parallel site calling.
5. VEP transcript-strand annotation.
6. Independent depth checks in all replicates.
7. Treated/control evidence matrix and optional WGS filtering.

### Caption

**Figure 1. Model 1 converts matched RNA-seq data into transcript-oriented RNA-editing evidence.** Paired-end reads from three editor-treated and three control libraries are aligned to GRCh38 and processed with GATK. REDItools2 identifies quality-supported substitutions, VEP supplies transcript orientation, and candidate coordinates are evaluated independently for sequencing depth in every replicate. The final evidence matrix retains strand-consistent C-to-U candidates that satisfy configurable replicate, depth, control, and genomic-variant criteria.

## Figure 2 — Strand-aware C-to-U definition

### Diagram text

```text
Positive-strand transcript
RNA:      C → U
Genome:   C → T

Negative-strand transcript
RNA:      C → U
Genome:   G → A
```

### Caption

**Figure 2. Genomic representation of transcript-level C-to-U editing.** On positive-strand transcripts, C-to-U editing appears as genomic C-to-T. On negative-strand transcripts, the reverse-complement representation is genomic G-to-A. The latter does not indicate biochemical editing of G.

## Figure 3 — Why independent depth confirmation is required

### Diagram concept

Show three control outcomes at the same coordinate:

```text
Control A: sufficient depth, no edited reads → informative negative
Control B: low depth                         → uncertain
Control C: no coverage                       → missing observation
```

### Caption

**Figure 3. A missing REDItools2 call is not automatically a negative observation.** Candidate-site depth is measured independently in every BAM so that an adequately sequenced negative sample can be distinguished from a low-coverage or uncovered site.

## Table 1 — Experimental design

| Condition | Replicate | Sample ID | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

## Table 2 — Main parameters

| Stage | Parameter | Value |
|---|---|---:|
| Reference | Genome assembly | GRCh38 primary assembly |
| Alignment | STAR mode | two-pass |
| REDItools2 | Edited-read threshold | 20 |
| Candidate depth | Minimum base quality | 30 |
| Candidate depth | Minimum mapping quality | 20 |
| Final filter | Minimum treated replicates | 3 of 3 |
| Final filter | Maximum control replicates | 0 of 3 |
| Final filter | Minimum depth | 20 in all six samples |
| Variant filtering | Matched WGS | optional but recommended |

## Table 3 — Evidence retained for every site

| Evidence field | Meaning |
|---|---|
| Chromosome and position | GRCh38 genomic coordinate |
| Reference and alternate base | Observed genomic substitution |
| VEP strand | Transcript orientation used for C-to-U interpretation |
| Treated called replicates | Number of treated libraries with a REDItools2 call |
| Control called replicates | Number of control libraries with a REDItools2 call |
| Candidate-site depth | Independent depth measurement in each BAM |
| Edited-read count | Reads supporting the alternate base |
| Editing frequency | Alternate-base fraction reported for the sample |
| WGS overlap | Whether the substitution overlaps an optional genomic variant set |

## Result figures to add later

The following panels require completed numerical outputs and are therefore not included yet:

1. six-sample UpSet plot of strand-consistent calls;
2. treated-versus-control candidate counts;
3. editing-frequency distribution;
4. depth distribution across retained candidates;
5. genomic or transcript-region annotation summary;
6. ranked candidate table with per-replicate support.
