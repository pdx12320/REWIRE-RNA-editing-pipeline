# Why we evaluated transcriptome-wide editing

Our REWIRE system recruits a cytidine deaminase to a selected RNA sequence. Reporter editing demonstrates that the designed editor is active, but it does not show whether the same construct also changes endogenous RNAs. We therefore built an RNA-seq workflow to identify transcript-oriented C-to-U candidate events across the transcriptome.

The analysis compares three editor-treated libraries with three matched controls. Replicate consistency is central to the design: a mismatch found in only one library may reflect sequencing error, alignment uncertainty, or stochastic low-level noise. A stronger candidate should be supported independently in treated replicates and evaluated against control coverage at the same coordinate.
