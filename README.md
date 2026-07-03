# REWIRE Model 1 — RNA-editing evidence pipeline

This repository contains the complete implementation of **Model 1**, the RNA-editing evidence pipeline used in the REWIRE iGEM project.

Model 1 converts six paired-end RNA-seq libraries into an auditable, transcript-oriented evidence matrix for candidate C-to-U editing sites.

## Main links

- **Copy-ready iGEM Wiki page:** [`wiki/Model1_RNA_editing_evidence_pipeline.md`](wiki/Model1_RNA_editing_evidence_pipeline.md)
- **Editable Wiki figures:** [`wiki/assets/`](wiki/assets/)
- **Sample manifest:** [`config/samples.tsv`](config/samples.tsv)
- **Executable scripts:** [`scripts/`](scripts/)
- **Technical notes:** [`docs/`](docs/)

The Wiki page links back to this README so readers can inspect the code and methodological details directly on GitHub.

---

## Scientific question

Reporter editing demonstrates that REWIRE can act at the designed target, but it does not establish transcriptome-wide specificity. Model 1 therefore asks:

1. Which substitutions are reproducibly detected in editor-treated samples?
2. Which substitutions are consistent with transcript-level C-to-U editing?
3. Which candidates remain after control, depth, strand, and optional genomic-variant filtering?

---

## Dataset

| Group | Replicate | Sample | SRA accession |
|---|---:|---|---|
| Treated | 1 | CU517_GC_T1 | SRR27885768 |
| Treated | 2 | CU517_GC_T2 | SRR27885766 |
| Treated | 3 | CU517_GC_T3 | SRR27885765 |
| Control | 1 | CU517_GC_C1 | SRR27885767 |
| Control | 2 | CU517_GC_C2 | SRR27885764 |
| Control | 3 | CU517_GC_C3 | SRR27885763 |

---

## Workflow

![Model 1 workflow](wiki/assets/model1_workflow.svg)

```text
SRA paired-end RNA-seq
→ FASTQ conversion
→ STAR two-pass alignment to GRCh38
→ read-group validation
→ GATK MarkDuplicates
→ GATK SplitNCigarReads
→ coverage-aware parallel REDItools2 calling
→ union substitution VCF
→ VEP transcript-strand annotation
→ transcript-oriented C-to-U interpretation
→ depth confirmation in all six samples
→ treated/control comparison
→ optional matched HEK293T WGS filtering
```

---

## Repository layout

```text
config/
  samples.tsv

scripts/
  download_sra_fastq.py
  run_star_alignment.py
  run_gatk_preprocessing.py
  generate_reditools_coverage_limited.sh
  run_reditools_all_samples.sh
  reditools_union_to_vcf.py
  run_vep_annotation.py
  build_candidate_depth_tables.sh
  filter_utils.py
  filter_calls.py
  filter_c_to_u_and_compare.py

wiki/
  Model1_RNA_editing_evidence_pipeline.md
  assets/
    model1_workflow.svg
    strand_orientation.svg
    evidence_logic.svg
    filtering_funnel.svg
    contig_fix.svg
    wetlab_drylab_loop.svg
    dbtl.svg

docs/
  installation, commands, troubleshooting, outputs, and limitations

results/
  placeholders only; final numerical results are not committed yet
```

---

## Quick start

Define the shared paths:

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
STAR_INDEX=/data/ydx/igem/STAR_index
REDITOOLS=/data/ydx/igem/REDItools2
MANIFEST=config/samples.tsv
```

### 1. Download SRA and convert to FASTQ

```bash
python3 scripts/download_sra_fastq.py \
  --project "$PROJECT" \
  --manifest "$MANIFEST" \
  --threads 16
```

### 2. STAR alignment

```bash
python3 scripts/run_star_alignment.py \
  --project "$PROJECT" \
  --star-index "$STAR_INDEX" \
  --manifest "$MANIFEST" \
  --threads 50
```

### 3. GATK RNA preprocessing

```bash
python3 scripts/run_gatk_preprocessing.py \
  --project "$PROJECT" \
  --reference "$REF" \
  --manifest "$MANIFEST" \
  --java-options=-Xmx16g
```

### 4. REDItools2

```bash
conda activate reditools2_py2

nohup bash scripts/run_reditools_all_samples.sh \
  "$PROJECT" \
  "$REF" \
  "$REDITOOLS" \
  "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools.log" 2>&1 &
```

The last three values are:

```text
30  MPI processes
8   simultaneous coverage jobs
8   compression threads
```

### 5. Union VCF and BED

```bash
python3 scripts/reditools_union_to_vcf.py \
  --manifest "$MANIFEST" \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"
```

### 6. VEP transcript-strand annotation

```bash
python3 scripts/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache /path/to/vep_cache
```

### 7. Candidate-site depth in all six samples

```bash
bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  "$MANIFEST"
```

### 8. Final evidence matrix

```bash
python3 scripts/filter_c_to_u_and_compare.py \
  --manifest "$MANIFEST" \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vep "$PROJECT/vep/CU5.17_EGFP_GC.vep.tsv" \
  --depth-dir "$PROJECT/candidate_depth" \
  --output-dir "$PROJECT/final" \
  --min-treated-reps 3 \
  --max-control-called-reps 0 \
  --min-depth-all-reps 20
```

Optional matched WGS filtering:

```text
--wgs-vcf /path/to/HEK293T.filtered.vcf.gz
```

---

## Key methodological decisions

### Strict REDItools2 discovery

```text
-S       output positions containing observed edits
-me 20   require at least 20 editing events at a reported position
-q 20    default minimum mapping quality
-bq 30   default minimum base quality
```

`-me 20` is an edited-read threshold, not a total-depth threshold.

### Transcript-oriented C-to-U interpretation

![Strand interpretation](wiki/assets/strand_orientation.svg)

```text
positive-strand transcript: genomic C→T
negative-strand transcript: genomic G→A
```

### Independent depth confirmation

A missing REDItools2 call is interpreted together with independently measured candidate-site depth in every BAM.

![Evidence logic](wiki/assets/evidence_logic.svg)

### Conservative default rule

![Filtering funnel](wiki/assets/filtering_funnel.svg)

```text
called in all three treated replicates
called in no control replicate
candidate-site depth ≥20 in all six samples
strand consistent with transcript-level C-to-U editing
absent from the optional matched WGS variant set
```

---

## REDItools2 contig-name fix

GRCh38 supplementary contigs include version suffixes such as `.1` and `.2`. The original temporary-file parser removed these suffixes and caused the final merge to fail.

![Contig parsing fix](wiki/assets/contig_fix.svg)

The corrected parser is:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

Only the final `.gz` extension is removed; the exact contig identifier is preserved.

---

## Quality-control checks

```bash
samtools quickcheck -v sample.bam
gzip -t sample.txt.gz
tabix -l sample.txt.gz
zcat sample.txt.gz | head
```

The workflow also checks read groups, per-sample coverage directories, interval completion, reference ordering, and the presence of tabix indexes.

---

## Wiki figures

All diagrams are stored as editable SVG files under `wiki/assets/`:

- complete workflow;
- transcript-strand interpretation;
- depth and replicate evidence logic;
- evidence filtering funnel;
- REDItools2 contig parsing fix;
- Wet Lab–Dry Lab feedback loop;
- Design–Build–Test–Learn.

They can be downloaded directly from GitHub and uploaded to the iGEM Wiki media system without redrawing.

---

## Results status

The methods, source code, and Wiki figures are included. Numerical result files, final site counts, and result-specific plots are intentionally excluded until all six samples complete the same workflow and pass the same integrity checks.

The original repository state remains available on branch `backup-before-paper-pipeline-20260704`.
