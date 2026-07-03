# Pipeline design

The project uses three treated RNA-seq replicates and three matched controls. Reads are aligned to GRCh38 with STAR. GATK read-group validation, duplicate marking, and SplitNCigarReads are applied before site discovery.

REDItools2 analyzes each sample independently. A per-position coverage map is generated first so the MPI program can divide the genome into balanced intervals. GL and KI names in the log are supplementary GRCh38 contigs and are expected.

Observed substitutions are combined into a union variant table and annotated with VEP transcript strand. Transcript-level C-to-U editing is represented as genomic C-to-T on positive-strand transcripts and genomic G-to-A on negative-strand transcripts.

Candidate positions are checked for sequencing depth in every replicate. The default treatment-specific definition requires reproducibility in all treated samples, no call in controls, adequate depth in all six samples, and no overlap with the optional matched HEK293T genomic variant set.

The final matrix retains per-sample call status, depth, edited-read count, and editing rate so that thresholds can be reviewed or changed.
