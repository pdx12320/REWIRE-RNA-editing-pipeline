# Model 1 — RNA-editing Evidence Pipeline

This folder contains the iGEM-ready documentation for **Model 1**, the RNA-editing evidence pipeline used in our REWIRE project.

Unlike a journal-style Methods section, the pages here follow an iGEM narrative:

```text
Why did we build it?
→ How did we design it?
→ What did we implement?
→ How did we test each step?
→ What did we learn?
→ How does it support the wet-lab project?
```

## Wiki-ready pages

1. [`Model1_RNA_Editing_Evidence_Pipeline.md`](Model1_RNA_Editing_Evidence_Pipeline.md)  
   Main Dry Lab page. This is the primary text to adapt into the team wiki.

2. [`Code_and_Reproducibility.md`](Code_and_Reproducibility.md)  
   Exact commands, script entry points, parameter explanations, and folder structure.

3. [`Figure_and_Layout_Plan.md`](Figure_and_Layout_Plan.md)  
   Suggested iGEM page layout, figures, cards, and result placeholders.

## What is included

- the biological question behind Model 1;
- the treated/control experimental design;
- the complete RNA-seq evidence workflow;
- code for every computational stage;
- quality-control checkpoints;
- the REDItools2 contig-name fix;
- strand-aware C-to-U interpretation;
- replicate and depth filtering logic;
- limitations and wet-lab validation plan.

## What is not included yet

Final result tables, numerical site counts, and completed result figures are intentionally excluded until all six samples finish the same workflow and pass the same quality checks.

## Main repository scripts

```text
scripts/download_sra_fastq.py
scripts/run_star_alignment.py
scripts/run_gatk_preprocessing.py
scripts/generate_reditools_coverage_limited.sh
scripts/run_reditools_all_samples.sh
scripts/reditools_union_to_vcf.py
scripts/run_vep_annotation.py
scripts/build_candidate_depth_tables.sh
scripts/filter_c_to_u_and_compare.py
```

The server run currently uses:

```text
/data/ydx/igem/run_cu517_egfp_gc_paper_pipeline_v2.sh
```

The modular scripts in the repository reproduce the same logic in a form that is easier to document and reuse.
