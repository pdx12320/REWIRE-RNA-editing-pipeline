# Results Summary

## Dataset

| Item | Value |
|---|---|
| Editor condition | CU5.15 |
| PUF target | AACGUCUAUA |
| Treated samples | SRR27885760, SRR27885759 |
| Reference genome | GRCh38 primary assembly |
| Annotation | GENCODE v50 primary assembly GTF |
| Minimum coverage | 20 |

## Candidate site counts

| Category | Count |
|---|---:|
| Control-corrected candidate sites | 33,983 |
| High-confidence candidate sites | 2,955 |
| Candidate sites with treated_edit_rate ≥ 0.10 | 5,995 |
| Candidate sites with treated_edit_rate ≥ 0.20 | 597 |
| Candidate sites with treated_edit_rate ≥ 0.50 | 6 |

## High-confidence candidate editing-rate distribution

| Treated edit-rate bin | Count |
|---|---:|
| 0.05–0.10 | 0 |
| 0.10–0.20 | 2,446 |
| 0.20–0.50 | 503 |
| ≥0.50 | 6 |

## Summary statistics for high-confidence candidates

| Metric | Value |
|---|---:|
| Median treated_edit_rate | 0.1404 |
| Mean treated_edit_rate | 0.1567 |
| Maximum treated_edit_rate | 0.9578 |
| Median corrected_edit_rate | 0.1400 |
| Maximum corrected_edit_rate | 0.9570 |
| Median treated_total_count | 77 |

## Main observation

The treated-control corrected candidate table contains 33,983 candidate C-to-U editing sites. Applying stricter high-confidence filters retained 2,955 sites.

Among the high-confidence candidates, most sites are in the 0.10–0.20 treated editing-rate range. A smaller number of sites show 0.20–0.50 editing, and only 6 sites show treated_edit_rate ≥ 0.50.

These results support using the dataset as a candidate editing classification or ranking dataset, with high-confidence candidate sites as positives and clean background C sites as negatives.
