#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def save_pipeline_overview(out):
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.axis('off')
    steps = [
        'SRA files',
        'FASTQ extraction',
        'STAR alignment',
        'Sorted BAM',
        'C-equivalent site scanning',
        'Base counting',
        'Replicate count merging',
        'Treated-control correction',
        'Candidate editing sites',
        'LAMAR labels',
    ]
    y0 = 0.92
    dy = 0.085
    for i, step in enumerate(steps):
        y = y0 - i * dy
        ax.text(0.5, y, step, ha='center', va='center', fontsize=13,
                bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='black'))
        if i < len(steps) - 1:
            ax.annotate('', xy=(0.5, y - 0.045), xytext=(0.5, y - 0.025),
                        arrowprops=dict(arrowstyle='->', lw=1.5))
    ax.set_title('RNA-seq to Control-corrected C-to-U Editing Labels', fontsize=16, pad=18)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def save_strand_definition(out):
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.axis('off')
    ax.set_title('Strand-aware C-to-U Editing Signal Definition', fontsize=16, pad=18)

    ax.text(0.05, 0.72, 'Positive-strand transcript', fontsize=14, weight='bold')
    ax.text(0.08, 0.57, 'transcript C = genomic C', fontsize=13)
    ax.text(0.08, 0.45, 'C-to-U editing appears as genomic C→T', fontsize=13)

    ax.text(0.55, 0.72, 'Negative-strand transcript', fontsize=14, weight='bold')
    ax.text(0.58, 0.57, 'transcript C = genomic G', fontsize=13)
    ax.text(0.58, 0.45, 'C-to-U editing appears as genomic G→A', fontsize=13)

    ax.text(0.5, 0.15,
            'G→A on the negative strand is the genomic-coordinate representation of transcript-level C→U editing.',
            ha='center', fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def save_treated_vs_control(candidates, out):
    df = pd.read_csv(candidates, sep='\t')
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(df['control_edit_rate'], df['treated_edit_rate'], s=6, alpha=0.35)
    ax.plot([0, 1], [0, 1], linestyle='--', linewidth=1)
    ax.set_xlabel('Control edit rate')
    ax.set_ylabel('Treated edit rate')
    ax.set_title('Treated vs Control Editing Rate')
    ax.set_xlim(0, min(1, max(0.25, df['control_edit_rate'].max() * 1.1)))
    ax.set_ylim(0, min(1, max(0.25, df['treated_edit_rate'].max() * 1.1)))
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def save_edit_rate_distribution(candidates, out):
    df = pd.read_csv(candidates, sep='\t')
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(df['treated_edit_rate'], bins=50)
    ax.set_xlabel('Treated edit rate')
    ax.set_ylabel('Number of candidate sites')
    ax.set_title('Candidate Editing-rate Distribution')
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def save_candidate_level_barplot(candidates, out):
    df = pd.read_csv(candidates, sep='\t')
    bins = [0.05, 0.10, 0.20, 0.50, 1.01]
    labels = ['0.05–0.10', '0.10–0.20', '0.20–0.50', '≥0.50']
    counts = pd.cut(df['treated_edit_rate'], bins=bins, labels=labels, right=False).value_counts().reindex(labels).fillna(0)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(counts.index.astype(str), counts.values)
    ax.set_xlabel('Treated edit-rate bin')
    ax.set_ylabel('Number of candidate sites')
    ax.set_title('Candidate Editing Level')
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', default='results/CU5.15_EGFP_CC.control_corrected.PUF_candidates.tsv.gz')
    ap.add_argument('--outdir', default='figures')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    save_pipeline_overview(outdir / 'pipeline_overview.png')
    save_strand_definition(outdir / 'strand_definition.png')

    candidate_path = Path(args.candidates)
    if candidate_path.exists():
        save_treated_vs_control(candidate_path, outdir / 'treated_vs_control_scatter.png')
        save_edit_rate_distribution(candidate_path, outdir / 'edit_rate_distribution.png')
        save_candidate_level_barplot(candidate_path, outdir / 'candidate_level_barplot.png')
    else:
        print(f'Candidate table not found: {candidate_path}')
        print('Only schematic figures were generated.')


if __name__ == '__main__':
    main()
