# Editing Site Extraction Results Summary

Fill this page after running the treated and control pipeline.

## Dataset

| Item | Value |
|---|---|
| Editor condition | CU5.15 |
| PUF target | AACGUCUAUA |
| Treated samples | SRR IDs here |
| Mock/control samples | SRR IDs here |
| Reference genome | GRCh38 primary assembly |
| Annotation | GENCODE v50 primary assembly GTF |
| Minimum site coverage | 20 |

## Transcriptome-wide site counts

| Category | Count |
|---|---:|
| Treated covered C-equivalent sites | xxx |
| Control covered C-equivalent sites | xxx |
| Control-corrected candidate sites | xxx |
| High-confidence candidate sites | xxx |
| Strong candidate sites, edit_rate >= 0.5 | xxx |

## Candidate editing-rate distribution

| Treated edit-rate bin | Count |
|---|---:|
| 0.05–0.10 | xxx |
| 0.10–0.20 | xxx |
| 0.20–0.50 | xxx |
| >=0.50 | xxx |

## Main observation

Most candidate sites are expected to show low-to-moderate editing levels, while high-efficiency off-target candidate sites should be rare. Therefore, this dataset is more suitable for candidate editing classification/ranking than for precise editing-efficiency regression.

Suggested wording for Wiki:

> We observed detectable low-to-moderate transcriptome-wide candidate editing signals after treated-control correction, whereas high-efficiency candidate off-target events were rare.

## Suggested figures

1. Pipeline overview diagram.
2. Strand-aware C-to-U signal definition.
3. Treated vs control editing-rate scatter plot.
4. Candidate editing-rate bin plot.
5. Optional PUF target similarity distribution.
