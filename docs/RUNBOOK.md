# Runbook

The workflow is divided into prepare, download, align, gatk, reditools, vep, filter, and all stages. Existing outputs are reused unless the force option is supplied.

For the current analysis, STAR and GATK are complete. Run the REDItools2 stage from the dedicated Python 2 environment, using 30 MPI processes and 8 concurrent coverage jobs. After all six sample tables finish, run VEP strand annotation and then final treated-control filtering.

The current project uses the GRCh38 primary assembly, a matching STAR index, and the REDItools2 repository installed on the analysis server. The default final criteria require a call in all three treated replicates, no call in control replicates, and at least 20 reads of candidate-site depth in every replicate. A matched HEK293T WGS VCF should be added when available.

Monitor coverage by checking active samtools depth processes and the size of the per-sample coverage directory. Monitor REDItools2 by following the progress lines in each sample log. Completion requires six compressed REDItools2 tables and six matching tabix index files.
