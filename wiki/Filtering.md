# Replicate and depth filtering

A site that is absent from a control call table may simply have insufficient sequencing coverage. We therefore measured candidate-site depth independently in all six BAM files. Our conservative default required a call in all three treated replicates, no call in any control replicate, and at least 20 reads of candidate-site depth in every sample. When matched HEK293T WGS data are available, overlapping genomic variants are removed before reporting the final treatment-specific table.
