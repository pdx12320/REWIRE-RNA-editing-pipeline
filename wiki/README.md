# HEK293T RNA-editing workflow

- [Workflow overview and flowchart](../README.md)
- [Complete filtering specification](../pipeline/CURRENT_HEK293T_WORKFLOW.md)
- [Recorded settings](../pipeline/config/current_hek293t.json)
- [Available tools](../pipeline/README.md)

Six editors: PUF10, PUF12, 132D, E72A, GVE and SNE. Four biological replicates per group; Control/mock, APOBEC-only and PUF-only supply background evidence. The final positive-only export requires ALT >=20 and depth >=100 in every treated replicate. Editing rate is the median of the four replicate rates. No 10% cutoff and no negative class are applied.
