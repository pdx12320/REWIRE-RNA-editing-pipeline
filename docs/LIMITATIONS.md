# Interpretation and limitations

A computational candidate is not automatically a validated biological off-target. RNA-seq mismatches can arise from genomic variants, alignment ambiguity, sequencing error, endogenous RNA modification, or true editor activity.

The REDItools2 setting requiring at least 20 edited reads is intentionally stringent. It favors well-covered, higher-frequency signals and can miss real low-frequency events. A later sensitivity analysis may use lower edited-read thresholds with stronger artifact controls.

Transcript strand is not always unique. A genomic coordinate may overlap transcripts on opposite strands. The filtering code excludes conflicting strand annotations rather than assigning one orientation arbitrarily.

Absence from a control call table is interpretable only when the control has sufficient depth. The workflow therefore measures candidate-site depth in all six replicates independently of whether REDItools2 reported a call.

Matched HEK293T WGS filtering is optional in the code but important for strong biological claims. Without it, the final output should be described as an RNA-derived candidate list rather than a definitive set of editing-only events.

High-priority sites should be confirmed using targeted amplicon sequencing, independent RNA-seq, Sanger sequencing for strong signals, or another orthogonal assay.
