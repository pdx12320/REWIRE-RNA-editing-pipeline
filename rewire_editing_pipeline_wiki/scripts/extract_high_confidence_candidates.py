#!/usr/bin/env python3
import argparse
import gzip

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, help='Control-corrected candidate table')
    ap.add_argument('--output', required=True, help='Output high-confidence candidate table')
    ap.add_argument('--min-treated-rate', type=float, default=0.10)
    ap.add_argument('--min-corrected-rate', type=float, default=0.10)
    ap.add_argument('--max-control-rate', type=float, default=0.02)
    ap.add_argument('--max-fisher-q', type=float, default=0.05)
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep='\t')

    mask = (
        (df['treated_edit_rate'] >= args.min_treated_rate) &
        (df['corrected_edit_rate'] >= args.min_corrected_rate) &
        (df['control_edit_rate'] <= args.max_control_rate)
    )

    if 'fisher_q' in df.columns:
        mask = mask & (df['fisher_q'] <= args.max_fisher_q)

    out = df.loc[mask].copy()
    out.to_csv(args.output, sep='\t', index=False, compression='gzip')

    print(f'input_rows={len(df)}')
    print(f'output_rows={len(out)}')
    print(f'output_file={args.output}')


if __name__ == '__main__':
    main()
