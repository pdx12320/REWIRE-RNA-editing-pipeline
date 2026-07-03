# Methods

## RNA-seq alignment and preprocessing

Paired-end RNA-seq reads from three editor-treated samples and three matched controls were aligned to the GRCh38 primary assembly using STAR in two-pass mode. Coordinate-sorted BAM files were generated with sample-level read-group metadata. BAM headers were checked before GATK processing; missing read groups were repaired when necessary. PCR duplicates were marked with GATK MarkDuplicates, and SplitNCigarReads was applied to process spliced RNA alignments before site-level analysis.

## REDItools2 candidate discovery

Each sample was analyzed independently with the MPI implementation of REDItools2. A per-position coverage map was generated before site calling so genomic intervals could be distributed across workers according to sequencing depth. Coverage was calculated for every sequence in the GRCh38 FASTA index, including supplementary GL and KI contigs. The number of concurrent `samtools depth` processes was limited to reduce disk contention on the shared server.

REDItools2 was run in strict mode with a minimum of 20 edited reads. Parallel interval files were sorted according to the GRCh38 FASTA index, merged, compressed with bgzip, and indexed with tabix. Contig version suffixes such as `.1` and `.2` were preserved during filename parsing.

## Strand-aware C-to-U interpretation

Substitutions observed across all six REDItools2 tables were combined into a union VCF and annotated with Ensembl VEP in offline GRCh38 mode. Transcript-level C-to-U editing was defined as genomic C-to-T on positive-strand transcripts or genomic G-to-A on negative-strand transcripts. Coordinates with contradictory transcript-strand annotations were treated as ambiguous and excluded.

## Replicate depth and treated-control comparison

Union candidate coordinates were queried independently in every SplitNCigarReads BAM using a minimum base quality of 30 and mapping quality of 20. This allowed a missing REDItools2 call to be distinguished from insufficient sequencing coverage.

The conservative default evidence rule required a candidate to be called in all three treated replicates, called in no control replicate, and covered by at least 20 reads in all six samples. When a matched HEK293T WGS VCF is available, overlapping genomic variants are removed. The final matrix retains per-sample call status, candidate-site depth, REDItools2 depth, edited-read count, and editing frequency.
