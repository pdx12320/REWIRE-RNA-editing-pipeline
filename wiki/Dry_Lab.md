# Dry Lab

This page describes the transcriptome-wide analysis developed for the CU5.17_EGFP_GC RNA-editing dataset. Three treated RNA-seq replicates are compared with three matched controls.

The workflow aligns paired-end reads to GRCh38, preprocesses RNA alignments, calls substitutions with REDItools2, annotates transcript strand with VEP, and compares replicate-level signals. Candidate-site depth is measured independently in every sample so that a missing call is not confused with missing sequencing coverage.

The analysis interprets transcript-level C-to-U events as genomic C-to-T on positive-strand transcripts and genomic G-to-A on negative-strand transcripts. The final candidate matrix preserves per-sample call status, depth, edited-read count, and editing frequency.

Detailed method and result-page text is provided in the other files in this folder.
