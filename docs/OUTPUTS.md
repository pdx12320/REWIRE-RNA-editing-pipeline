# Output files

The external project directory is organized into raw data, aligned BAM files, GATK outputs, REDItools2 files, VEP annotations, candidate-depth tables, final calls, and logs.

## Per-sample REDItools2 tables

Six compressed tables are expected under `reditools/tables/`, one for each treated or control replicate. Each table also has a tabix index. The columns report chromosome, position, reference base, inferred strand, quality-filtered depth, mean quality, A/C/G/T counts, observed substitution, and editing frequency. The genomic-DNA fields remain empty when no DNA BAM is supplied directly to REDItools2.

## Candidate-depth tables

The `candidate_depth/` directory contains one table per sample. Each row stores chromosome, position, and quality-filtered depth at a union candidate site. These files distinguish a true negative call from an unsequenced site.

## Annotation files

The `vcf/` directory contains the union substitution VCF and candidate-position BED. The `vep/` directory contains offline VEP annotations used to determine transcript orientation.

## Final result tables

The `final/` directory contains per-sample transcript-oriented C-to-U calls, a complete six-sample site matrix, a treated-consensus table, and a treatment-specific table.

The principal final file is `CU5.17_EGFP_GC.treatment_specific.tsv.gz`. It records coordinates, reference and alternate bases, VEP strand, treated and control replicate counts, all-replicate depth status, optional WGS overlap, and per-replicate measurements.

No numerical result tables are committed until all six samples have completed the same workflow and passed integrity checks.
