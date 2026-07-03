# Model 1 — Suggested iGEM Page Layout

This layout turns the pipeline into a project story rather than a manuscript-style Methods page.

## 1. Hero

**Model 1 — RNA-editing Evidence Pipeline**  
*From six RNA-seq libraries to replicate-supported evidence for REWIRE activity*

Use one simple visual:

```text
RNA-seq → evidence pipeline → ranked candidates → validation
```

## 2. The problem

Use two short columns:

**What we knew**

- REWIRE edited the designed reporter.
- RNA-seq data were available for treated and control samples.

**What we still needed to know**

- whether editing signals appeared elsewhere;
- whether those signals were reproducible;
- whether controls had enough coverage;
- which sites should be validated first.

## 3. Our design

Show six sample cards:

```text
T1  T2  T3     C1  C2  C3
```

Explain that replicates test reproducibility and controls provide background context.

## 4. Full workflow figure

```text
SRA
→ STAR
→ GATK
→ REDItools2
→ VEP strand annotation
→ all-sample depth check
→ treated/control evidence matrix
```

Keep the biological purpose larger than the software names.

## 5. Pipeline cards

Use five cards:

1. **Map** — align RNA-seq reads to GRCh38.
2. **Prepare** — process read groups, duplicates, and splice-aware alignments.
3. **Call** — identify supported substitutions with REDItools2.
4. **Orient** — convert genomic calls into transcript-level C-to-U evidence.
5. **Compare** — evaluate every site across all six samples.

Each card should show:

- one-sentence purpose;
- one key output;
- one quality-control check;
- a link to the full code.

## 6. Strand figure

```text
Positive transcript: genomic C→T
Negative transcript: genomic G→A
```

Both represent transcript-level C-to-U editing.

## 7. Evidence logic figure

Use three example cards:

```text
No control call + high control depth = informative negative
No control call + zero control depth = missing observation
Call in only one treated replicate = not reproducible
```

## 8. Engineering story

Show the GL/KI contig parsing problem and the fix:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

Present this as a debugging contribution, not as a hidden implementation note.

## 9. Filtering funnel

```text
all substitutions
→ strand-consistent C-to-U
→ 3/3 treated support
→ 0/3 control calls
→ depth ≥20 in all samples
→ optional WGS exclusion
→ validation candidates
```

Leave numerical counts blank until the complete run is finished.

## 10. Wet Lab–Dry Lab loop

```text
RNA-seq
→ Model 1 evidence
→ candidate ranking
→ targeted validation
→ improved labels for future models
```

## 11. Design–Build–Test–Learn

Use four equal cards:

- **Design:** six-sample evidence rules;
- **Build:** modular RNA-seq pipeline;
- **Test:** file, depth, strand, and merge checks;
- **Learn:** absence needs coverage, strand changes representation, exact contig names matter.

## 12. Results placeholder

Reserve sections for:

- run completion and QC;
- replicate overlap;
- treated/control comparison;
- candidate frequency and depth;
- ranked candidates;
- experimental validation.

Do not upload result tables or provisional counts before all six samples complete the same analysis.
