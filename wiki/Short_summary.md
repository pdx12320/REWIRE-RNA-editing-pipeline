# Model 1 summary for the iGEM Wiki

To evaluate the transcriptome-wide specificity of our REWIRE cytidine editor, we developed a reproducible RNA-seq evidence pipeline using three treated and three matched control libraries. Reads are aligned to GRCh38 with STAR, processed with GATK, and analyzed with parallel REDItools2. Candidate substitutions are then annotated with VEP so that transcript-level C-to-U events can be interpreted correctly as genomic C-to-T on positive-strand transcripts or genomic G-to-A on negative-strand transcripts.

A missing call is not automatically treated as evidence of absence. Candidate-site depth is measured independently in all six samples, and the final evidence matrix retains replicate-level call status, depth, edited-read support, and editing frequency. Our conservative default requires support in all treated replicates, no call in controls, sufficient depth in every sample, and removal of matched genomic variants when WGS data are available.

Model 1 provides the evidence-generation layer for the REWIRE dry lab. It converts raw sequencing data into auditable computational RNA-editing candidates that can support downstream biological interpretation and predictive modeling.
