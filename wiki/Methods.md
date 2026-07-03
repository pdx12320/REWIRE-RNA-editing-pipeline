# Methods

RNA-seq reads were aligned to GRCh38 with STAR and processed with GATK before REDItools2 site discovery. VEP transcript-strand annotation was then used to interpret genomic substitutions in transcript orientation.

Candidate positions were checked for sequencing depth in every treated and control replicate. The final table stores per-sample calls, depth, edited-read support, and editing frequency.
