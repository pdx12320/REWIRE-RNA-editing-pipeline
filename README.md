# REWIRE RNA-editing pipeline

Current workflow: transfected human HEK293T RNA-seq, adapter trimming, GRCh38 plus GTF alignment, quality-filtered editing evidence, three-control background filtering, and six-editor comparison. This is a **positive-only computational screening workflow**.

## Current study

Six treatments: PUF10, PUF12, 132D, E72A, GVE and SNE; four biological replicates each. Controls: mock/Control, APOBEC-only and PUF-only; four replicates each. Total: 36 libraries. Raw sample identifiers may retain `GVD`; report labels use **GVE**.

- BQ >30 and MAPQ >30 mean integer thresholds **31**.
- Initial positive screen: four treated replicates each ALT >=5, AF >=0.005; coverage >=20 in every sample of the declared analysis universe.
- Final exported sites and off-target figures: **each of the four treated replicates ALT >=20 and depth >=100**.
- Editing rate: **median of four replicate ALT/(REF+ALT) rates**, not pooled ALT divided by pooled depth.
- No 10% editing-rate or 10-percentage-point increase cutoff.
- No negative class is generated.

## Workflow

```mermaid
flowchart TD
    raw["Paired FASTQ and manifest"] --> trim["Adapter trimming and QC"]
    trim --> align["Human GRCh38 plus GTF"]
    align --> bam["STAR and GATK preprocessing"]
    bam --> call["REDItools candidate calling"]
    call --> recount["BQ31 MAPQ31 NH1 recount"]
    recount --> universe{"Analysis scope"}
    universe -->|"PUF12 priority"| sixteen["16 sample coverage"]
    universe -->|"Six editor comparison"| thirtySix["36 sample coverage"]
    sixteen --> filter["SNP and background filters"]
    thirtySix --> filter
    filter --> positive["Four replicate positive screen"]
    positive --> finalSites["Each replicate ALT20 depth100"]
    finalSites --> export["Tables and separate figures"]
```

See [complete workflow and screening rules](pipeline/CURRENT_HEK293T_WORKFLOW.md), [machine-readable settings](pipeline/config/current_hek293t.json), and [pipeline entry point](pipeline/README.md).

## Export the final subset

The following **postprocessing command does not perform alignment, calling or background filtering**. Its inputs must already be background/SNP-filtered positives and their quality-filtered replicate counts.

```bash
python3 pipeline/scripts/rna/export_alt20_depth100.py \
  --positive PUF12_both_positive.tsv \
  --counts PUF12_positive_replicate_counts.tsv \
  --samples PUF12_1 PUF12_2 PUF12_3 PUF12_4 \
  --output PUF12_positive_ALT20_depth100.tsv
```

Run the same exporter for each treatment. It preserves all input site columns and adds per-replicate depth, ALT count and rate. It rejects duplicate/incomplete evidence and existing output files.

## Add strand to final positives

```bash
python3 pipeline/scripts/rna/export_positive_strand.py \
  --positive PUF12_positive_ALT20_depth100.tsv \
  --output PUF12_positive_with_strand.tsv
```

This annotation-only wrapper accepts already transcript-strand-compatible C-to-U positives: genomic `C>T` maps to `strand=+`, and `G>A` maps to `strand=-`. It does not infer gene strand from arbitrary variants or read orientation. Upstream transcript annotation must already have established compatibility.

It preserves coordinates, genomic alleles, editing rates, all other columns, row order and already oriented sequences. If `sequence_101nt` is present, it must contain 101 unambiguous DNA bases with C at index 50; it is **not reverse-complemented again**. Unsupported alleles, duplicate sites, conflicting existing strands and existing output files are rejected. No background filtering or threshold changes are performed.

Run the standard-library tests with `python3 -m unittest discover -s tests -v`.

## Interpretation and status

The independent PUF12 release uses a 16-sample callable universe (four PUF12 plus twelve controls). The final six-group comparison uses a 36-sample common callable universe. Their counts need not match. Do not combine these as if their denominators were identical.

The deployed final six-group rerun, background filtering and figure generation completed on **2026-09-24**. This records computational completion, not biological validation or a portable end-to-end reproduction test. C388 is an apparent T/(C+T) measurement confounded by endogenous reference T. Computational positive sites are not experimentally proven off-target events. Result tables and figures are not included in this code update.

Raw reads, BAMs, private sample paths, credentials and machine-specific deployment wrappers are not distributed here.
