# Current HEK293T positive-only workflow

Specification date: 2026-09-22. This supersedes the historical CU5.17 workflow as the repository entry point. Existing historical scripts/data retain their original semantics.

## Inputs and reference

Use paired-end FASTQ, an explicit sample/group/replicate manifest, human GRCh38 primary-assembly FASTA and matching GENCODE v50 GTF, STAR index, offline GRCh38 VEP cache, and four SNP catalogues: 293T_RTG, 293T_DMBR_CG, 293T_14_CG and 293T_CG. HEK293T is human, not mouse. Verify GTF contig names and coordinate bounds against the FASTA.

There are four biological replicates per group. Groups are Control/mock, APOBEC-only, PUF-only, PUF10, PUF12, 132D, E72A, GVE and SNE. The input alias GVD maps to display name GVE without renaming immutable source files.

## Preprocessing

1. Inspect dependencies, FASTQ existence/checksums, disk, CPU and memory. Log missing files and errors; never invent results or silently skip missing samples.
2. Trim paired-end Illumina adapters with fastp. Current adapters: R1 `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`; R2 `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`. Preserve JSON/HTML QC. Disable additional quality filtering, length filtering and poly-G trimming in this specific run; quality selection occurs during read counting.
3. STAR two-pass Basic, human GRCh38, matching GTF via `--sjdbGTFfile`, `--sjdbOverhang 149` for 150-nt reads, coordinate-sorted BAM, `NH HI AS nM MD`, unique MAPQ60 and explicit read groups.
4. GATK MarkDuplicates then SplitNCigarReads. Audit and restore uniqueness tags when the latter drops them: reconstruction from MAPQ60 is allowed only after verifying source NH1/MAPQ60 equivalence. Never assume this equivalence for other aligner settings.
5. REDItools2 candidate discovery uses `-S -me 5 -bq 31 -q 31 -men 1`. **-me is total non-reference support, not target-ALT support.** Keep -men1 because -men5 can reject an entire column when any other substitution has fewer than five reads. Target-ALT >=5 is enforced explicitly in the recount-based screen.
6. Union calls, VEP annotation and direct candidate recount in every included BAM. Use BQ>=31, MAPQ>=31, NH=1; exclude unmapped, duplicate, secondary, supplementary, QC-failed, deleted/skipped and non-ACGT observations. Retain the pileup overlap behavior and software version in provenance; do not mix count implementations.

This study evaluates APOBEC C-to-U: genomic C>T and reverse-complement G>A, normalized to a 101-nt C-centered sequence (zero-based center 50). A-to-I is not included in these final APOBEC tables. Canonical nuclear chromosomes are chr1–chr22, chrX and chrY.

## Quality, SNP and sequence gates

Depth in the reported allele rate is **REF count + target ALT count**, not total A/C/G/T coverage. Require depth>=20 in every sample of the declared universe; missing evidence is not zero editing.

For the four provided SNP files, preserve the verified source assembly (hg18/build36), select PASS biallelic SNVs with non-reference genotype, lift to GRCh38, validate REF and exclude the union by exact CHROM/POS/REF/ALT. This is not exclusion of all variants at the same position and not a claim of exhaustive donor-specific WGS filtering.

Require a complete unambiguous 101-nt context centered on C, Shannon entropy>=1.2, maximum homopolymer run<20 and dinucleotide-repeat fraction<0.8. This is a low-complexity screen, not an external mappability assay. Retain gene/transcript annotations separately; multiple sites in one gene remain distinct sites.

## Background screen

Use all three groups: Control/mock, APOBEC-only and PUF-only. In each group require:

- median of four replicate rates <=0.02;
- zero of four formal caller detections;
- median absolute deviation <=0.02.

Also require APOBEC-only minus Control median <=0.02 and PUF-only minus Control median <=0.02. These are low-background gates, **not a requirement for zero ALT reads**.

For each treatment, perform a one-sided pooled-read Fisher test against each control group; take the maximum of the three p-values per site, then BH-adjust over the eligible callable universe. Require FDR<=0.05. Pooled-read FDR is a screening statistic, not a biological-replicate hypothesis test.

## Positive and final export gates

| Gate | Initial positive | Final off-target export/figure |
|---|---|---|
| Treated formal detection | 4/4 replicates | retained |
| Treated target ALT reads | >=5 in each replicate | **>=20 in each replicate** |
| Treated allelic fraction | >=0.005 in each replicate | retained |
| Treated depth | >=20 in each replicate | **>=100 in each replicate** |
| Treated rate MAD | <=0.05 | retained |
| SNP, sequence, background and FDR | pass | retained |

No threshold on editing rate >10% or treatment-background increase >10 percentage points is applied. A high depth does not imply high ALT support.

Per-replicate rate is ALT/(REF+ALT). Site `editing_rate` is the median of four rates; with four replicates this is the mean of the two middle values. It is **not** sum(ALT)/sum(depth). The exported editing rate is raw, not numerically background-subtracted.

## Scope and execution order

- Priority: complete PUF12 and all twelve controls, call/recount/filter over this **16-sample universe**, export the PUF12 result.
- Continue remaining treatments, then repeat union/recount/filter using **all 36 samples** as the common coverage universe.
- Apply final ALT20/depth100 gates to each treatment's four replicates and generate six-group figures. Controls retain their depth20/background criteria; ALT20 is not imposed on controls.
- Preserve the independent PUF12 output. Changing 16 to 36 samples can change both callability and the FDR testing universe, so the final PUF12 subset can differ.

Use isolated output directories, timestamped logs, command and environment records, lock files and per-stage validated checkpoints. Reuse only completed preprocessing with unchanged parameters. Changed calling/recount thresholds require new outputs, not reuse of old threshold markers. Retain incomplete outputs for inspection.

The deployed resource envelope was two concurrent samples, 24 STAR threads each, 12 MPI ranks per calling sample and eight recount workers. These are resource settings, not universal recommendations; validate memory and disk before choosing concurrency.

## Targets, figures and limitations

Report C295, C388 and C871 independently of off-target positivity. For target plots require depth>=100 in all four treatment replicates; do not impose ALT20 on targets because it would hide low target efficiency. Missing/insufficient-depth values are NA, not zero.

Current GRCh38 target coordinates: chr19:44908591 (C295), chr19:44908684 (C388), chr19:44909167 (C871). C388 has endogenous reference T, so T/(C+T) is **apparent efficiency** and cannot separate edited transgene from endogenous T. Construct-aware evidence is needed for unconfounded efficiency. Exclude C388 from off-target counts; C295/C871 are currently not excluded and must be disclosed.

Output separate editing-rate distribution (site medians plus C388 marker), C295/C388/C871 comparison, motif logos/composition and chromosome distribution. Use GVE labels, no abcd labels or unnecessary explanatory text inside plots; put caveats in captions. Motifs are descriptive frequencies, not matched-background enrichment. Preserve source tables, random jitter seed and visual review status.

## What this update distributes

This repository update supplies the current specification, flowchart, settings and a portable final-subset exporter. It does **not** turn historical CU5.17 runners into a newly validated portable 36-sample launcher. Deployment-specific orchestration remains outside the repository until its paths/dependencies are parameterized and integration-tested. No raw sequencing data, new negative set, or final six-group count claim is published.
