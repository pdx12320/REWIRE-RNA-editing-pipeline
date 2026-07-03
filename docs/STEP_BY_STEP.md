# Step-by-step execution

1. Download the six SRA accessions and convert each library to paired FASTQ files.
2. Align paired reads to the GRCh38 STAR index in two-pass mode and write sample read groups.
3. Run GATK read-group repair when needed, followed by MarkDuplicates and SplitNCigarReads.
4. Activate the REDItools2 Python 2 environment and run `run_reditools_all_samples.sh` with the project path, reference FASTA, REDItools2 directory, Python interpreter, 30 MPI processes, 8 coverage jobs, and 8 compression threads.
5. Run `reditools_union_to_vcf.py` to create the union substitution VCF and candidate BED.
6. Annotate the union VCF with offline VEP for GRCh38, retaining Uploaded_variation, Location, Allele, and STRAND.
7. Run `build_candidate_depth_tables.sh` to measure candidate-site depth in all six SplitNCigarReads BAM files.
8. Run `filter_c_to_u_and_compare.py` with the manifest, REDItools2 tables, VEP table, candidate-depth directory, and optional HEK293T WGS VCF.
9. Copy final tables and figures into the repository results area only after integrity checks.

The operational settings used on the current server are 50 alignment threads, 30 MPI processes, 8 simultaneous coverage jobs, base quality 30, mapping quality 20, and minimum all-replicate depth 20.
