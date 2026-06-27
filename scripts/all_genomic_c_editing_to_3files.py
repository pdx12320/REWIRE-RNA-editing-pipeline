#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import pysam


def parse_gtf_exons(gtf_path):
    exons = defaultdict(list)
    with open(gtf_path) as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'exon':
                continue
            chrom, _, _, start, end, _, strand, _, attrs = parts
            start0 = int(start) - 1
            end0 = int(end)
            gene_id = gene_name = transcript_id = ''
            for item in attrs.split(';'):
                item = item.strip()
                if not item or ' ' not in item:
                    continue
                key, val = item.split(' ', 1)
                val = val.strip().strip('"')
                if key == 'gene_id':
                    gene_id = val
                elif key == 'gene_name':
                    gene_name = val
                elif key == 'transcript_id':
                    transcript_id = val
            exons[(chrom, strand)].append((start0, end0, gene_id, gene_name, transcript_id))

    merged = defaultdict(list)
    for key, items in exons.items():
        items = sorted(items, key=lambda x: (x[0], x[1]))
        cur_s = cur_e = None
        payload = []
        for s, e, gid, gname, tid in items:
            if cur_s is None:
                cur_s, cur_e = s, e
                payload = [(s, e, gid, gname, tid)]
            elif s <= cur_e:
                cur_e = max(cur_e, e)
                payload.append((s, e, gid, gname, tid))
            else:
                merged[key].append((cur_s, cur_e, payload))
                cur_s, cur_e = s, e
                payload = [(s, e, gid, gname, tid)]
        if cur_s is not None:
            merged[key].append((cur_s, cur_e, payload))
    return merged


def open_bams(bam_dir, srrs):
    bams = []
    for srr in srrs:
        p1 = Path(bam_dir) / f'{srr}.split.bam'
        p2 = Path(bam_dir) / f'{srr}.Aligned.sortedByCoord.out.bam'
        if p1.exists():
            path = p1
        elif p2.exists():
            path = p2
        else:
            raise FileNotFoundError(f'Missing BAM for {srr}: {p1} or {p2}')
        bams.append(pysam.AlignmentFile(str(path), 'rb'))
    return bams


def count_site(bams, chrom, pos0, edited_base, min_baseq):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    for bam in bams:
        try:
            pileup = bam.pileup(
                chrom,
                pos0,
                pos0 + 1,
                truncate=True,
                stepper='nofilter',
                min_base_quality=min_baseq,
                ignore_overlaps=False,
                ignore_orphans=False,
            )
        except ValueError:
            continue
        for col in pileup:
            if col.reference_pos != pos0:
                continue
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                base = pr.alignment.query_sequence[pr.query_position].upper()
                if base not in counts:
                    base = 'N'
                counts[base] += 1
    total = sum(counts.values())
    edited = counts.get(edited_base, 0)
    rate = edited / total if total else 0.0
    return counts, edited, total, rate


def site_annotation(payload, pos0):
    gene_ids, gene_names, transcript_ids = set(), set(), set()
    for s, e, gid, gname, tid in payload:
        if s <= pos0 < e:
            if gid:
                gene_ids.add(gid)
            if gname:
                gene_names.add(gname)
            if tid:
                transcript_ids.add(tid)
    return ';'.join(sorted(gene_names)), ';'.join(sorted(gene_ids)), ';'.join(sorted(transcript_ids))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--srr', required=True, help='Comma-separated SRR IDs')
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--editor', default='')
    ap.add_argument('--puf-target', default='')
    ap.add_argument('--ref', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--bam-dir', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--min-coverage', type=int, default=20)
    ap.add_argument('--min-baseq', type=int, default=20)
    ap.add_argument('--edited-rate-cutoff', type=float, default=0.05)
    ap.add_argument('--strong-rate-cutoff', type=float, default=0.5)
    args = ap.parse_args()

    srrs = [x.strip() for x in args.srr.split(',') if x.strip()]
    ref = pysam.FastaFile(args.ref)
    bams = open_bams(args.bam_dir, srrs)
    merged_exons = parse_gtf_exons(args.gtf)

    all_path = f'{args.out_prefix}.all_C_sites.genomic.cov{args.min_coverage}.tsv.gz'
    edited_path = f'{args.out_prefix}.edited_sites.genomic.cov{args.min_coverage}.rate005.tsv.gz'
    strong_path = f'{args.out_prefix}.strong_sites.genomic.cov{args.min_coverage}.rate05.tsv.gz'

    header = [
        'sample_id', 'editor', 'puf_target_seq',
        'chrom', 'strand', 'genomic_pos_1based',
        'transcript_base', 'genomic_ref_base', 'edited_base_in_genome',
        'edited_count', 'total_count', 'edit_rate',
        'A_count', 'C_count', 'G_count', 'T_count', 'N_count',
        'gene_names', 'gene_ids', 'transcript_ids',
    ]

    n_intervals = sum(len(v) for v in merged_exons.values())
    scanned = all_n = edited_n = strong_n = 0

    with gzip.open(all_path, 'wt') as fall, gzip.open(edited_path, 'wt') as fedit, gzip.open(strong_path, 'wt') as fstrong:
        for fh in (fall, fedit, fstrong):
            fh.write('\t'.join(header) + '\n')

        for (chrom, strand), intervals in merged_exons.items():
            for start, end, payload in intervals:
                scanned += 1
                if scanned % 500 == 0:
                    print(f'scanned_intervals={scanned}/{n_intervals} all={all_n} edited>={args.edited_rate_cutoff}={edited_n} strong>={args.strong_rate_cutoff}={strong_n}', flush=True)

                try:
                    seq = ref.fetch(chrom, start, end).upper()
                except Exception:
                    continue

                for i, base in enumerate(seq):
                    pos0 = start + i
                    if strand == '+':
                        if base != 'C':
                            continue
                        genomic_ref_base = 'C'
                        edited_base = 'T'
                    else:
                        if base != 'G':
                            continue
                        genomic_ref_base = 'G'
                        edited_base = 'A'

                    counts, edited, total, rate = count_site(bams, chrom, pos0, edited_base, args.min_baseq)
                    if total < args.min_coverage:
                        continue

                    gene_names, gene_ids, transcript_ids = site_annotation(payload, pos0)
                    row = [
                        args.sample_id, args.editor, args.puf_target,
                        chrom, strand, str(pos0 + 1),
                        'C', genomic_ref_base, edited_base,
                        str(edited), str(total), f'{rate:.6g}',
                        str(counts['A']), str(counts['C']), str(counts['G']), str(counts['T']), str(counts['N']),
                        gene_names, gene_ids, transcript_ids,
                    ]
                    line = '\t'.join(row) + '\n'
                    fall.write(line)
                    all_n += 1
                    if rate >= args.edited_rate_cutoff:
                        fedit.write(line)
                        edited_n += 1
                    if rate >= args.strong_rate_cutoff:
                        fstrong.write(line)
                        strong_n += 1

    for bam in bams:
        bam.close()
    ref.close()
    print(f'Done. all={all_n} edited={edited_n} strong={strong_n}', flush=True)


if __name__ == '__main__':
    main()
