#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import sqlite3
import tempfile


def open_text(path, mode='rt'):
    return gzip.open(path, mode, newline='') if str(path).endswith('.gz') else open(path, mode, newline='')


def to_int(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default


def to_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def rc(seq):
    tab = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(tab)[::-1].upper()


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def best_puf_match(seq, puf):
    seq = seq.upper()
    puf = puf.upper().replace('U', 'T')
    k = len(puf)
    best_d, best_pos, best_kmer = None, None, ''
    if len(seq) < k:
        return '', '', ''
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i:i + k]
        if 'N' in kmer:
            continue
        d = hamming(kmer, puf)
        if best_d is None or d < best_d:
            best_d, best_pos, best_kmer = d, i, kmer
    if best_d is None:
        return '', '', ''
    return str(best_d), str(best_pos), best_kmer


def motif_annotation(ref, chrom, strand, pos1, puf, flank, max_mismatch):
    pos0 = int(pos1) - 1
    start = max(0, pos0 - flank)
    end = pos0 + flank + 1
    try:
        seq = ref.fetch(chrom, start, end).upper()
    except Exception:
        return {'best_puf_hamming': '', 'best_puf_rel_start': '', 'best_puf_kmer': '', 'pass_puf_like': ''}

    center = pos0 - start
    if strand == '-':
        seq = rc(seq)
        center = len(seq) - 1 - center

    d, best_pos, kmer = best_puf_match(seq, puf)
    if d == '':
        return {'best_puf_hamming': '', 'best_puf_rel_start': '', 'best_puf_kmer': '', 'pass_puf_like': ''}

    rel = int(best_pos) - center
    pass_like = str(int(d) <= max_mismatch)
    return {'best_puf_hamming': d, 'best_puf_rel_start': str(rel), 'best_puf_kmer': kmer, 'pass_puf_like': pass_like}


def build_control_db(control_path, db_path, chunksize=200000):
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute('PRAGMA journal_mode=OFF')
    cur.execute('PRAGMA synchronous=OFF')
    cur.execute('PRAGMA temp_store=MEMORY')
    cur.execute('''
        CREATE TABLE control (
            chrom TEXT,
            strand TEXT,
            pos INTEGER,
            edited INTEGER,
            total INTEGER,
            rate REAL,
            PRIMARY KEY (chrom, strand, pos)
        )
    ''')

    n = 0
    batch = []
    with open_text(control_path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for r in reader:
            batch.append((
                r['chrom'],
                r['strand'],
                to_int(r['genomic_pos_1based']),
                to_int(r['edited_count']),
                to_int(r['total_count']),
                to_float(r['edit_rate']),
            ))
            n += 1
            if len(batch) >= chunksize:
                cur.executemany('INSERT OR REPLACE INTO control VALUES (?,?,?,?,?,?)', batch)
                con.commit()
                batch.clear()
                print(f'loaded_control_rows={n}', flush=True)
        if batch:
            cur.executemany('INSERT OR REPLACE INTO control VALUES (?,?,?,?,?,?)', batch)
            con.commit()
            print(f'loaded_control_rows={n}', flush=True)

    cur.execute('CREATE INDEX IF NOT EXISTS idx_control_key ON control(chrom, strand, pos)')
    con.commit()
    return con


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--treated', required=True)
    ap.add_argument('--control', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-treated-coverage', type=int, default=20)
    ap.add_argument('--min-control-coverage', type=int, default=20)
    ap.add_argument('--min-treated-rate', type=float, default=0.05)
    ap.add_argument('--max-control-rate', type=float, default=0.02)
    ap.add_argument('--min-delta', type=float, default=0.05)
    ap.add_argument('--min-corrected-rate', type=float, default=0.05)
    ap.add_argument('--tmp-db', default='')
    ap.add_argument('--ref', default='')
    ap.add_argument('--puf-target', default='')
    ap.add_argument('--flank', type=int, default=50)
    ap.add_argument('--max-puf-mismatch', type=int, default=10)
    args = ap.parse_args()

    ref = None
    if args.ref and args.puf_target:
        import pysam
        ref = pysam.FastaFile(args.ref)

    if args.tmp_db:
        db_path = args.tmp_db
        if os.path.exists(db_path):
            os.remove(db_path)
    else:
        tmp = tempfile.NamedTemporaryFile(prefix='control_editing_', suffix='.sqlite', delete=False)
        db_path = tmp.name
        tmp.close()

    print('Building sqlite index for control:', args.control, flush=True)
    con = build_control_db(args.control, db_path)
    cur = con.cursor()

    base_fields = [
        'sample_id', 'editor', 'puf_target_seq',
        'chrom', 'strand', 'genomic_pos_1based',
        'transcript_base', 'genomic_ref_base', 'edited_base_in_genome',
        'treated_edited_count', 'treated_total_count', 'treated_edit_rate',
        'control_edited_count', 'control_total_count', 'control_edit_rate',
        'delta_edit_rate', 'corrected_edit_rate',
        'treated_unedited_count', 'control_unedited_count',
        'A_count', 'C_count', 'G_count', 'T_count', 'N_count',
        'gene_names', 'gene_ids', 'transcript_ids',
    ]
    motif_fields = ['best_puf_hamming', 'best_puf_rel_start', 'best_puf_kmer', 'pass_puf_like'] if ref is not None else []

    print('Scanning treated and writing:', args.out, flush=True)
    n_read = n_merged = n_kept = 0
    with open_text(args.treated, 'rt') as fin, open_text(args.out, 'wt') as fout:
        reader = csv.DictReader(fin, delimiter='\t')
        writer = csv.DictWriter(fout, fieldnames=base_fields + motif_fields, delimiter='\t')
        writer.writeheader()

        for r in reader:
            n_read += 1
            chrom = r['chrom']
            strand = r['strand']
            pos = to_int(r['genomic_pos_1based'])

            cur.execute('SELECT edited,total,rate FROM control WHERE chrom=? AND strand=? AND pos=?', (chrom, strand, pos))
            hit = cur.fetchone()
            if hit is None:
                if n_read % 500000 == 0:
                    print(f'read={n_read} merged={n_merged} candidates={n_kept}', flush=True)
                continue

            n_merged += 1
            c_edited, c_total, c_rate = int(hit[0]), int(hit[1]), float(hit[2])
            t_edited = to_int(r['edited_count'])
            t_total = to_int(r['total_count'])
            t_rate = to_float(r['edit_rate'])

            if t_total < args.min_treated_coverage or c_total < args.min_control_coverage:
                continue

            delta = t_rate - c_rate
            corrected = 0.0
            if c_rate < 1:
                corrected = max(0.0, min(1.0, (t_rate - c_rate) / (1.0 - c_rate)))

            if t_rate < args.min_treated_rate:
                continue
            if c_rate > args.max_control_rate:
                continue
            if delta < args.min_delta:
                continue
            if corrected < args.min_corrected_rate:
                continue

            out = {
                'sample_id': r.get('sample_id', ''),
                'editor': r.get('editor', ''),
                'puf_target_seq': r.get('puf_target_seq', ''),
                'chrom': chrom,
                'strand': strand,
                'genomic_pos_1based': pos,
                'transcript_base': r.get('transcript_base', 'C'),
                'genomic_ref_base': r.get('genomic_ref_base', ''),
                'edited_base_in_genome': r.get('edited_base_in_genome', ''),
                'treated_edited_count': t_edited,
                'treated_total_count': t_total,
                'treated_edit_rate': f'{t_rate:.6g}',
                'control_edited_count': c_edited,
                'control_total_count': c_total,
                'control_edit_rate': f'{c_rate:.6g}',
                'delta_edit_rate': f'{delta:.6g}',
                'corrected_edit_rate': f'{corrected:.6g}',
                'treated_unedited_count': t_total - t_edited,
                'control_unedited_count': c_total - c_edited,
                'A_count': r.get('A_count', ''),
                'C_count': r.get('C_count', ''),
                'G_count': r.get('G_count', ''),
                'T_count': r.get('T_count', ''),
                'N_count': r.get('N_count', ''),
                'gene_names': r.get('gene_names', ''),
                'gene_ids': r.get('gene_ids', ''),
                'transcript_ids': r.get('transcript_ids', ''),
            }

            if ref is not None:
                out.update(motif_annotation(ref, chrom, strand, pos, args.puf_target, args.flank, args.max_puf_mismatch))

            writer.writerow(out)
            n_kept += 1
            if n_read % 500000 == 0:
                print(f'read={n_read} merged={n_merged} candidates={n_kept}', flush=True)

    if ref is not None:
        ref.close()
    con.close()
    if not args.tmp_db:
        try:
            os.remove(db_path)
        except Exception:
            pass

    print('Done.', flush=True)
    print(f'read={n_read} merged={n_merged} candidates={n_kept}', flush=True)


if __name__ == '__main__':
    main()
