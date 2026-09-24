"""Add transcript/editing strand to already strand-compatible C-to-U positives.

This is an annotation-only export, not a caller or a strand inference tool for
arbitrary variants. Inputs must already have passed transcript-aware annotation.
Genomic alleles and pre-oriented sequence contexts are preserved verbatim.
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

KEY = ('chromosome', 'position', 'ref', 'alt')
STRANDS = {('C', 'T'): '+', ('G', 'A'): '-'}


def annotate(fields, rows):
    if not fields or len(fields) != len(set(fields)):
        raise ValueError('Missing or duplicate column names')
    if not set(KEY).issubset(fields):
        raise ValueError('Required columns: ' + ', '.join(KEY))
    output_fields = list(fields)
    if 'strand' not in fields:
        output_fields.insert(output_fields.index('position') + 1, 'strand')
    result, seen = [], set()
    for line, original in enumerate(rows, start=2):
        if set(original) != set(fields) or any(v is None for v in original.values()):
            raise ValueError(f'Line {line}: malformed TSV row')
        if not original['chromosome'] or not original['position'].isdigit() or int(original['position']) < 1:
            raise ValueError(f'Line {line}: expected a nonempty chromosome and positive 1-based position')
        key = tuple(original[k] for k in KEY)
        if key in seen:
            raise ValueError(f'Line {line}: duplicate site {key}')
        seen.add(key)
        strand = STRANDS.get((original['ref'], original['alt']))
        if strand is None:
            raise ValueError(f'Line {line}: expected genomic C>T or G>A, got {key}')
        if 'strand' in fields and original['strand'] != strand:
            raise ValueError(f'Line {line}: existing strand conflicts with the C-to-U convention')
        if 'sequence_101nt' in fields:
            seq = original['sequence_101nt']
            if len(seq) != 101 or seq[50] != 'C' or set(seq) - set('ACGT'):
                raise ValueError(f'Line {line}: expected an already oriented unambiguous 101-nt C-centered sequence')
        result.append({**original, 'strand': strand})
    return output_fields, result


def export(source, destination):
    with Path(source).open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields, rows = annotate(reader.fieldnames, list(reader))
    # Validate the entire input before creating any output. Never replace files.
    with Path(destination).open('x', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    return Counter(row['strand'] for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--positive', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    counts = export(args.positive, args.output)
    print(f"Exported {sum(counts.values())} sites: +={counts['+']}, -={counts['-']}")


if __name__ == '__main__':
    main()
