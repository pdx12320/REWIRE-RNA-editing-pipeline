# Module 1 Wiki

- [Results and methods](01_Module1_Wiki_EN.md)
- [Design–Build–Test–Learn](02_Module1_DBTL_EN.md)
- [Complete filtering specification](../pipeline/CURRENT_HEK293T_WORKFLOW.md)
- [Recorded settings](../pipeline/config/current_hek293t.json)
- [Available tools](../pipeline/README.md)

The two English pages present the six-editor comparison for Wiki readers. Keep this directory intact when importing it into Obsidian or a Wiki build so that relative image and data links resolve.

## Figures and supporting data

The candidate-editing section uses the original [editing distribution](assets/03_editing_distribution.png), including site-level points, site medians, candidate counts and apparent C388 rates. It does not substitute a results table or separate count and median-summary charts.

Other panels show the [illustrated workflow](assets/01_module1_workflow_illustrated.png), [target comparison](assets/04_target_comparison.png), [sequence logos](assets/05_motif_logos.png) and [chromosome profiles](assets/06_chromosome_distribution.png).

[Displayed summaries](assets/data/displayed_summary.tsv), [screening summaries](assets/data/six_group_summary.tsv), [chromosome normalization data](assets/data/chromosome_reference_normalized.tsv) and [reference-intersection counts](assets/data/sample_scope_audit.tsv) accompany the figures. These are summary data, not raw sequencing or per-site distribution data.

## Editable workflow source

The workflow is available as [PowerPoint](figure-sources/illustrated-workflow/figure.pptx) and [PDF](figure-sources/illustrated-workflow/figure.pdf). Text, boxes and connectors are native PowerPoint objects. Three independently generated transparent illustrations are embedded as bitmap components; their [prompts](figure-sources/illustrated-workflow/prompts.json) and PNG assets are retained alongside the presentation. Molecular illustrations are conceptual and do not depict measured sequences.

## Interpretation

Candidate sites are computational screening outputs, not experimentally proven off-target events. C388 is an apparent signal affected by endogenous reference T. It is excluded from candidate counts, whereas C295 and C871 are not systematically excluded. An empty E72A candidate set does not establish an active, highly specific editor. Analysis thresholds and statistical methods are unchanged by this documentation update.
