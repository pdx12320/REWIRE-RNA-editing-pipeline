# Installation

This workflow requires STAR, samtools, GATK 4.4.0.0, REDItools2, Open MPI, bgzip, tabix, VEP, Python 3, and a separate Python 2.7 environment for REDItools2.

Use GRCh38 consistently for the FASTA, STAR index, GATK dictionary, and VEP cache. The REDItools2 environment must provide mpi4py, pysam, sortedcontainers, psutil, and netifaces. Compile mpi4py against the same MPI implementation used by mpirun.

Pass the exact Python 2 interpreter to the pipeline with the reditools-python option. An optional matched HEK293T WGS VCF can be supplied during final filtering.
