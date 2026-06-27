#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def save_pipeline_overview(out):
    fig, ax = plt.subplots(figsize=(10, 7.5))
    ax.axis('off')
    ax.set_title('RNA-seq to C-to-U Candidate Editing Sites', fontsize=16, pad=20)
    steps = [
        'SRA files',
        'FASTQ extraction',
        'STAR alignment to GRCh38',
        'Sorted BAM files',
        'GENCODE exon C-equivalent site scanning',
        'Base counting at each site',
        'Replicate count merging',
        'Treated-control correction',
        'Candidate C-to-U editing sites',
        'High-confidence candidate subset',
    ]
    y0 = 0.92
    dy = 0.085
    for i, step in enumerate(steps):
        y = y0 - i * dy
        ax.text(0.5, y, step, ha='center', va='center', fontsize=11,
                bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='black'))
        if i < len(steps) - 1:
            ax.annotate('', xy=(0.5, y - 0.055), xytext=(0.5, y - 0.025),
                        arrowprops=dict(arrowstyle='->', lw=1.2))
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def save_strand_definition(out):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis('off')
    ax.set_title('Strand-aware C-to-U Signal Definition', fontsize=16, pad=20)

    ax.text(0.08, 0.75, 'Positive-strand transcript', fontsize=13, weight='bold')
    ax.text(0.10, 0.60, 'Transcript C = genomic C', fontsize=12)
    ax.text(0.10, 0.48, 'C-to-U signal = genomic C→T', fontsize=12)

    ax.text(0.56, 0.75, 'Negative-strand transcript', fontsize=13, weight='bold')
    ax.text(0.58, 0.60, 'Transcript C = genomic G', fontsize=12)
    ax.text(0.58, 0.48, 'C-to-U signal = genomic G→A', fontsize=12)

    ax.text(0.5, 0.18,
            'G→A is the genomic-coordinate representation of transcript-level C→U on the negative strand.',
            ha='center', fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def save_treated_vs_control(df, out):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(df['control_edit_rate'], df['treated_edit_rate'], s=8, alpha=0.35)
    upper = min(1.0, max(0.25, df['treated_edit_rate'].max() * 1.05))
    ax.plot([0, upper], [0, upper], linestyle='--', linewidth=1)
    ax.set_xlabel('Control edit rate')
    ax.set_ylabel('Treated edit rate')
    ax.set_title('High-confidence Candidates: Treated vs Control')
    ax.set_xlim(0, max(0.03, df['control_edit_rate'].max() * 1.2))
    ax.set_ylim(0, upper)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def save_edit_rate_distribution(df, out):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(df['treated_edit_rate'], bins=50)
    ax.set_xlabel('Treated edit rate')
    ax.set_ylabel('Number of high-confidence candidate sites')
    ax.set_title('Editing-rate Distribution')
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def save_candidate_level_barplot(df, out):
    bins = [0.05, 0.10, 0.20, 0.50, 1.01]
    labels = ['0.05–0.10', '0.10–0.20', '0.20–0.50', '≥0.50']
    counts = pd.cut(df['treated_edit_rate'], bins=bins, labels=labels, right=False).value_counts().reindex(labels).fillna(0).astype(int)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(counts.index.astype(str), counts.values)
    ax.set_xlabel('Treated edit-rate bin')
    ax.set_ylabel('Number of high-confidence candidate sites')
    ax.set_title('Candidate Editing Level')
    for i, v in enumerate(counts.values):
        ax.text(i, v + max(counts.values) * 0.02, str(v), ha='center', va='bottom', fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', default='results/CU5.15_EGFP_CC.high_confidence_candidates.tsv.gz')
    ap.add_argument('--outdir', default='figures')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    save_pipeline_overview(outdir / 'pipeline_overview.png')
    save_strand_definition(outdir / 'strand_definition.png')

    candidate_path = Path(args.candidates)
    if not candidate_path.exists():
        print(f'Candidate table not found: {candidate_path}')
        print('Only schematic figures were generated.')
        return

    df = pd.read_csv(candidate_path, sep='\t')
    save_treated_vs_control(df, outdir / 'treated_vs_control_scatter.png')
    save_edit_rate_distribution(df, outdir / 'edit_rate_distribution.png')
    save_candidate_level_barplot(df, outdir / 'candidate_level_barplot.png')


if __name__ == '__main__':
    main()
