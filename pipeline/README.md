# Pipeline guide

See [current HEK293T workflow](CURRENT_HEK293T_WORKFLOW.md) for the full preprocessing, calling and screening specification and [settings](config/current_hek293t.json) for recorded parameters.

## Distributed tools

- [Final subset exporter](scripts/rna/export_alt20_depth100.py): export already background/SNP-filtered positives, requiring ALT >=20 and depth >=100 in each of four replicates. See the [usage example](../README.md#export-the-final-subset).
- [Candidate base counter](scripts/rna/pileup_candidate_base_counts.py): general recount utility. Its defaults are not the current protocol; explicitly provide MAPQ31/BQ31 and audit all other read/count semantics before use. This utility alone does not perform the three-control screen.
- [Counting environment](env/rna_counting.yml), [catalogue environment](env/genomic_catalogue.yml), and [REDItools Python environment](env/reditools2_py2.yml). These are component environments, not a complete deployment lockfile.

The repository distributes specifications and postprocessing utilities, not a portable end-to-end launcher. Machine-specific orchestration and raw data are not included. No negative-label or model-training workflow is distributed.
