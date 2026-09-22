"""Final subset export from ALREADY background/SNP-filtered positives."""
import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

KEY = ('chromosome', 'position', 'ref', 'alt')

def read(path):
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return reader.fieldnames, list(reader)

def select(sites, counts, samples):
    if len(samples) != 4 or len(set(samples)) != 4:
        raise ValueError('Exactly four distinct biological replicates required')
    keys = [tuple(r[k] for k in KEY) for r in sites]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate positive site')
    evidence = defaultdict(dict)
    for row in counts:
        key = tuple(row[k] for k in KEY)
        sample = row['sample']
        if sample not in samples:
            raise ValueError('Unexpected sample: ' + sample)
        if sample in evidence[key]:
            raise ValueError('Duplicate site/sample count')
        depth, alt = int(row['depth']), int(row['alt_reads'])
        if not 0 <= alt <= depth or depth <= 0:
            raise ValueError('Invalid allele counts')
        evidence[key][sample] = (depth, alt)
    selected = []
    for site, key in zip(sites, keys):
        records = evidence[key]
        if set(records) != set(samples):
            raise ValueError('Missing replicate evidence: ' + str(key))
        if not all(d >= 100 and a >= 20 and a / d >= .005 for d, a in records.values()):
            continue
        row = dict(site)
        rates = [records[s][1] / records[s][0] for s in samples]
        rate = float(site['editing_rate'])
        if not math.isfinite(rate) or abs(rate - statistics.median(rates)) > 1e-8:
            raise ValueError('Site rate is not median of four replicate ALT/depth rates')
        for sample in samples:
            d, a = records[sample]
            row.update({sample + '_depth': d, sample + '_alt_reads': a,
                        sample + '_editing_rate': a / d})
        selected.append(row)
    return selected

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--positive', required=True)
    p.add_argument('--counts', required=True)
    p.add_argument('--samples', nargs=4, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    fields, sites = read(args.positive)
    _, counts = read(args.counts)
    selected = select(sites, counts, args.samples)
    fields = fields + [s + '_' + k for s in args.samples for k in ('depth', 'alt_reads', 'editing_rate')]
    with args.output.open('x', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(selected)
    print(f'{len(sites)} background-filtered positives -> {len(selected)} final sites')

if __name__ == '__main__':
    main()
