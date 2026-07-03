# Model 1 — iGEM Wiki Layout and Figure Plan

This file translates the Model 1 content into a visual iGEM page. The page should feel like a project story rather than a manuscript.

---

## Recommended page order

### 1. Hero section

**Title:** Model 1 — RNA-editing Evidence Pipeline  
**Subtitle:** From six RNA-seq libraries to replicate-supported evidence for REWIRE activity

Suggested hero visual:

```text
RNA reads → genome alignment → editing evidence → wet-lab validation
```

Keep the hero short. Do not place software-version tables or long parameter lists above the first figure.

---

### 2. The problem

Use a two-column block.

**Left: What the wet lab can measure**

- editing at the designed reporter;
- expression of the REWIRE construct;
- selected validation sites.

**Right: What remained unknown**

- transcriptome-wide candidate activity;
- reproducibility across replicates;
- whether a missing control call reflects low coverage;
- which sites should be validated first.

Suggested callout:

> Reporter editing shows that REWIRE can work. Model 1 asks where else a reproducible editing signal may appear.

---

### 3. Our solution

## Figure 1 — Full workflow

Recommended horizontal stages:

```text
6 RNA-seq libraries
→ STAR
→ GATK
→ REDItools2
→ VEP strand annotation
→ six-sample evidence matrix
→ ranked candidates
```

Visual notes:

- treated samples use one visual family;
- control samples use another;
- merge both streams only at the replicate-comparison stage;
- place a small shield/check icon at quality-control checkpoints;
- show the wet-lab validation arrow returning from the final candidates.

Caption:

> Model 1 converts treated and control RNA-seq reads into transcript-oriented C-to-U candidate evidence. Each retained site is evaluated across all six samples rather than judged from one call table.

---

### 4. Experimental design card

Use six small sample cards rather than a dense table on the final webpage.

```text
Treated: T1, T2, T3
Control: C1, C2, C3
```

Under the cards, explain:

- replicates test reproducibility;
- controls estimate background;
- independent depth checks determine whether a missing call is informative.

---

### 5. Pipeline story

Use one illustrated section for each major stage rather than presenting every command in the main narrative.

Suggested cards:

1. **Map** — STAR aligns reads to GRCh38.
2. **Prepare** — GATK handles read groups, duplicates, and spliced alignments.
3. **Call** — REDItools2 finds supported RNA/reference mismatches.
4. **Orient** — VEP converts genomic substitutions into transcript-oriented C-to-U evidence.
5. **Compare** — all six replicates are combined into one evidence matrix.

Each card should contain:

- one sentence explaining the purpose;
- one key output;
- one quality-control check;
- a “View code” link to the repository.

Put the full commands in a collapsible code section or on a separate Software/Code page.

---

### 6. Strand explanation

## Figure 2 — Why C-to-U has two genomic representations

Two panels:

```text
Positive-strand transcript
RNA: C → U
Genome representation: C → T
```

```text
Negative-strand transcript
RNA: C → U
Genome representation: G → A
```

Caption:

> Genomic G-to-A on a negative-strand transcript is the reverse-complement representation of transcript-level C-to-U editing; it is not interpreted as biochemical editing of G.

This figure should appear before the final filtering rules.

---

### 7. Evidence logic

## Figure 3 — Why a missing call is not always a negative

Use three example boxes:

```text
Treated: called, depth 120
Control: not called, depth 130
Interpretation: informative treated-specific evidence
```

```text
Treated: called, depth 120
Control: not called, depth 0
Interpretation: control is uninformative
```

```text
Treated: called in only 1 of 3 replicates
Interpretation: not reproducible under the default rule
```

This is one of the strongest conceptual parts of Model 1 and should be shown visually.

---

### 8. Filtering rule

Use a funnel diagram rather than a paragraph.

```text
All REDItools2 substitutions
        ↓
Transcript-strand-consistent C-to-U
        ↓
Called in all three treated replicates
        ↓
Called in no control replicate
        ↓
Depth ≥20 in all six samples
        ↓
Optional matched-WGS exclusion
        ↓
Candidates for validation
```

Do not insert numerical site counts until the full analysis is complete.

---

### 9. Engineering and debugging

## Figure 4 — Contig parsing fix

Show one before/after example:

```text
Before
GL000009.2#100#500.gz
→ GL000009
→ reference mismatch and merge failure
```

```text
After
GL000009.2#100#500.gz
→ GL000009.2 | 100 | 500
→ correct reference ordering
```

Suggested code block:

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

Frame this as an engineering lesson: a reference-name mismatch can invalidate the final merge even after the expensive computation has finished.

---

### 10. Wet Lab–Dry Lab connection

## Figure 5 — Feedback loop

```text
Wet-lab RNA-seq
→ Model 1 evidence
→ candidate ranking
→ targeted validation
→ improved evidence and future model labels
```

Suggested text:

> Model 1 does not replace experimental validation. It reduces the search space and records why each site was prioritized.

---

### 11. Results section placeholder

Keep this section in the page structure, but do not upload result files yet.

Suggested future subsections:

- Run completion and QC;
- per-sample mismatch calls;
- strand-consistent C-to-U calls;
- treated replicate overlap;
- treated/control comparison;
- final candidate ranking;
- wet-lab validation status.

Suggested future figures:

1. six-sample UpSet plot;
2. editing-frequency distribution;
3. treated/control site matrix heat map;
4. depth distribution across all replicates;
5. ranked candidate table;
6. IGV snapshots for selected sites.

---

### 12. Limitations and next steps

Use two side-by-side cards.

**Current limitations**

- RNA-derived mismatches may include genomic variants or alignment artifacts;
- strict edited-read support can miss low-frequency activity;
- overlapping transcripts can create strand ambiguity;
- final sites require independent validation.

**Next steps**

- matched HEK293T WGS filtering;
- sensitivity analysis using alternative thresholds;
- targeted sequencing of high-priority sites;
- use validated evidence as labels for downstream sequence models.

---

## Recommended page components

Use:

- a short hero;
- one full-width workflow figure;
- alternating text/figure sections;
- colored treated and control sample cards;
- collapsible code blocks;
- small “Why this matters” callouts;
- a visible Design–Build–Test–Learn section;
- links to GitHub for complete scripts.

Avoid:

- a manuscript-style abstract;
- long uninterrupted Methods paragraphs;
- software citations mixed into the project story;
- raw server paths as the main visual focus;
- presenting unvalidated candidates as confirmed off-targets;
- adding provisional result numbers before the six-sample run is complete.
