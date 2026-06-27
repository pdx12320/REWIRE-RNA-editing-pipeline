#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import pysam


def parse_attrs(attr_text):
    out = {}
    for item in attr_text.split(';'):
        item = item.strip()
        if not item or ' ' not in item:
            continue
        key, val = item.split(' ', 1)
        out[key] = val.strip().strip('"')
    return out


def load_merged_exons(gtf):
    raw = defaultdict(list)
    with open(gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'exon':
                continue
            chrom, start, end, strand = p[0], int(p[3]) - 1, int(p[4]), p[6]
            a = parse_attrs(p[8])
            raw[(chrom, strand)].append((start, end, a.get('gene_id', ''), a.get('gene_name', ''), a.get('transcript_id', '')))

    merged = defaultdict(list)
    for key, items in raw.items():
        items.sort(key=lambda x: (x[0], x[1]))
        cur_s = cur_e = None
        payload = []
        for s, e, gid, gname, tid in items:
            if cur_s is None:
                cur_s, cur_e, payload = s, e, [(s, e, gid, gname, tid)]
            elif s <= cur_e:
                cur_e = max(cur_e, e)
                payload.append((s, e, gid, gname, tid))
            else:
                merged[key].append((cur_s, cur_e, payload))
                cur_s, cur_e, payload = s, e, [(s, e, gid, gname, tid)]
        if cur_s is not None:
            merged[key].append((cur_s, cur_e, payload))
    return merged


def open_bams(bam_dir, srrs):
    bams = []
    for srr in srrs:
        p1 = Path(bam_dir) / f'{srr}.split.bam'
        p2 = Path(bam_dir) / f'{srr}.Aligned.sortedByCoord.out.bam'
        if p1.exists():
            bams.append(pysam.AlignmentFile(str(p1), 'rb'))
        elif p2.exists():
            bams.append(pysam.AlignmentFile(str(p2), 'rb'))
        else:
            raise FileNotFoundError(f'Missing BAM for {srr}')
    return bams


def count_bases(bams, chrom, pos0, edited_base, min_baseq):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    for bam in bams:
        try:
            cols = bam.pileup(chrom, pos0, pos0 + 1, truncate=True, stepper='nofilter', min_base_quality=min_baseq)
        except ValueError:
            continue
        for col in cols:
            if col.reference_pos != pos0:
                continue
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                b = pr.alignment.query_sequence[pr.query_position].upper()
                counts[b if b in counts else 'N'] += 1
    total = sum(counts.values())
    edited = counts[edited_base]
    rate = edited / total if total else 0.0
    return counts, edited, total, rate


def anno_for_site(payload, pos0):
    gids, gnames, tids = set(), set(), set()
    for s, e, gid, gname, tid in payload:
        if s <= pos0 < e:
            if gid:
                gids.add(gid)
            if gname:
                gnames.add(gname)
            if tid:
                tids.add(tid)
    return ';'.join(sorted(gnames)), ';'.join(sorted(gids)), ';'.join(sorted(tids))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--srr', required=True)
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
    exons = load_merged_exons(args.gtf)

    all_file = f'{args.out_prefix}.all_C_sites.genomic.cov{args.min_coverage}.tsv.gz'
    edited_file = f'{args.out_prefix}.edited_sites.genomic.cov{args.min_coverage}.rate005.tsv.gz'
    strong_file = f'{args.out_prefix}.strong_sites.genomic.cov{args.min_coverage}.rate05.tsv.gz'

    header = ['sample_id','editor','puf_target_seq','chrom','strand','genomic_pos_1based','transcript_base','genomic_ref_base','edited_base_in_genome','edited_count','total_count','edit_rate','A_count','C_count','G_count','T_count','N_count','gene_names','gene_ids','transcript_ids']
    total_intervals = sum(len(v) for v in exons.values())
    scanned = n_all = n_edit = n_strong = 0

    with gzip.open(all_file, 'wt') as fa, gzip.open(edited_file, 'wt') as fe, gzip.open(strong_file, 'wt') as fs:
        for f in [fa, fe, fs]:
            f.write('\t'.join(header) + '\n')
        for (chrom, strand), intervals in exons.items():
            for start, end, payload in intervals:
                scanned += 1
                if scanned % 500 == 0:
                    print(f'scanned_intervals={scanned}/{total_intervals} all={n_all} edited>={args.edited_rate_cutoff}={n_edit} strong>={args.strong_rate_cutoff}={n_strong}', flush=True)
                try:
                    seq = ref.fetch(chrom, start, end).upper()
                except Exception:
                    continue
                for i, base in enumerate(seq):
                    pos0 = start + i
                    if strand == '+':
                        if base != 'C':
                            continue
                        ref_base, edit_base = 'C', 'T'
                    else:
                        if base != 'G':
                            continue
                        ref_base, edit_base = 'G', 'A'
                    counts, edited, cov, rate = count_bases(bams, chrom, pos0, edit_base, args.min_baseq)
                    if cov < args.min_coverage:
                        continue
                    gene_names, gene_ids, transcript_ids = anno_for_site(payload, pos0)
                    row = [args.sample_id,args.editor,args.puf_target,chrom,strand,str(pos0+1),'C',ref_base,edit_base,str(edited),str(cov),f'{rate:.6g}',str(counts['A']),str(counts['C']),str(counts['G']),str(counts['T']),str(counts['N']),gene_names,gene_ids,transcript_ids]
                    line = '\t'.join(row) + '\n'
                    fa.write(line)
                    n_all += 1
                    if rate >= args.edited_rate_cutoff:
                        fe.write(line)
                        n_edit += 1
                    if rate >= args.strong_rate_cutoff:
                        fs.write(line)
                        n_strong += 1
    for bam in bams:
        bam.close()
    ref.close()


if __name__ == '__main__':
    main()
