# Model 1 — RNA-editing Evidence Pipeline

Reporter editing tells us that REWIRE can work at the designed target. Model 1 asks the next question: **where else does a reproducible C-to-U signal appear?**

We compared three editor-treated RNA-seq libraries with three matched controls. Reads were aligned to GRCh38, prepared for RNA-aware site analysis, and scanned with parallel REDItools2. We then used transcript-strand annotation to interpret genomic C-to-T and G-to-A substitutions as transcript-level C-to-U evidence.

A missing call was not treated automatically as a negative result. Candidate-site depth was measured independently in all six samples, allowing us to distinguish a well-covered control without an edit call from a control that was never informative at that position.

The final evidence matrix retains replicate-level call status, sequencing depth, edited-read support, and editing frequency. Model 1 therefore acts as the evidence layer connecting wet-lab RNA-seq, candidate prioritization, experimental validation, and future sequence modeling.
