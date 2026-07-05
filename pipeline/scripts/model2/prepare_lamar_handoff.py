#!/usr/bin/env python3
"""Prepare a Lamar-ready metadata table from final REWIRE candidates.

The script summarizes treated-replicate editing rates and can optionally extract
a fixed-length, transcript-oriented sequence window from a GRCh38 FASTA.
"""

import argparse
import csv
import gzip
import re
import statistics
from pathlib import Path


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def as_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if text in {"", "NA", "NaN", "nan", "."}:
        return None
    return float(text)


def as_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def find_treated_rate_columns(fieldnames):
    cols = [
        name for name in fieldnames
        if name.endswith("_edit_rate") and re.search(r"(?:^|_)T[123](?:_|$)", name)
    ]
    if len(cols) != 3:
        expected = [f"CU517_GC_T{i}_edit_rate" for i in (1, 2, 3)]
        if all(name in fieldnames for name in expected):
            cols = expected
    if len(cols) != 3:
        raise SystemExit(
            "ERROR: expected three treated edit-rate columns; found: "
            + ", ".join(cols)
        )
    return sorted(cols)


def companion_column(rate_col, suffix, fieldnames):
    candidate = rate_col[:-len("_edit_rate")] + suffix
    return candidate if candidate in fieldnames else None


def load_fasta(reference):
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit(
            "ERROR: --reference requires pysam. Install it or omit --reference."
        ) from exc
    return pysam.FastaFile(str(reference))


def fetch_context(fasta, chrom, pos_1based, strand, flank):
    if fasta is None:
        return "", "not_extracted"

    names = set(fasta.references)
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append("chr" + chrom)
    contig = next((name for name in candidates if name in names), None)
    if contig is None:
        raise ValueError(f"contig not found in FASTA: {chrom}")

    start = pos_1based - flank - 1
    end = pos_1based + flank
    if start < 0:
        raise ValueError(f"window extends before contig start: {chrom}:{pos_1based}")

    seq = fasta.fetch(contig, start, end).upper()
    expected_len = 2 * flank + 1
    if len(seq) != expected_len:
        raise ValueError(
            f"window length {len(seq)} != {expected_len}: {chrom}:{pos_1based}"
        )
    if strand == -1:
        seq = reverse_complement(seq)
    if seq[flank] != "C":
        raise ValueError(
            f"transcript-oriented center is {seq[flank]}, not C: "
            f"{chrom}:{pos_1based}, strand={strand}"
        )
    return seq, "transcript_5to3_centered_on_edited_C"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Final treatment-specific TSV[.gz]")
    parser.add_argument("--output", required=True, help="Lamar handoff TSV[.gz]")
    parser.add_argument(
        "--reference",
        default="",
        help="Optional indexed GRCh38 FASTA for sequence extraction",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=50,
        help="Bases on each side of the edited C; default gives 101 nt",
    )
    args = parser.parse_args()

    fasta = load_fasta(args.reference) if args.reference else None

    with open_text(args.input) as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit("ERROR: input table has no header")
        fieldnames = list(reader.fieldnames)
        required = {"chrom", "position", "ref", "alt", "vep_strand"}
        missing = sorted(required - set(fieldnames))
        if missing:
            raise SystemExit("ERROR: missing columns: " + ", ".join(missing))

        rate_cols = find_treated_rate_columns(fieldnames)
        output_fields = [
            "site_id", "chrom", "position_grch38_1based",
            "genomic_ref", "genomic_alt", "transcript_strand",
            "transcript_ref", "transcript_alt",
            "treated_called_reps", "control_called_reps",
            "genomic_catalogue_overlap",
        ]
        for i, rate_col in enumerate(rate_cols, start=1):
            output_fields += [
                f"T{i}_coverage", f"T{i}_alt_count", f"T{i}_edit_rate"
            ]
        output_fields += [
            "median_treated_edit_rate", "mean_treated_edit_rate",
            "sd_treated_edit_rate", "min_treated_edit_rate",
            "max_treated_edit_rate", "range_treated_edit_rate",
            "label_recommendation", "label_caveat",
            "sequence_context", "sequence_orientation",
        ]

        rows = []
        for row in reader:
            chrom = row["chrom"]
            pos = as_int(row["position"])
            ref = row["ref"].upper()
            alt = row["alt"].upper()
            strand = as_int(row["vep_strand"])
            rates = [as_float(row[name]) for name in rate_cols]
            if any(value is None for value in rates):
                raise SystemExit(
                    f"ERROR: missing treated edit rate at {chrom}:{pos}"
                )

            sequence, orientation = fetch_context(
                fasta, chrom, pos, strand, args.flank
            )
            record = {
                "site_id": f"{chrom}:{pos}:{ref}>{alt}",
                "chrom": chrom,
                "position_grch38_1based": pos,
                "genomic_ref": ref,
                "genomic_alt": alt,
                "transcript_strand": strand,
                "transcript_ref": "C",
                "transcript_alt": "T",
                "treated_called_reps": as_int(row.get("treated_called_reps"), 3),
                "control_called_reps": as_int(row.get("control_called_reps"), 0),
                "genomic_catalogue_overlap": as_int(
                    row.get("genomic_catalogue_overlap", row.get("wgs_variant", 0))
                ),
            }
            for i, rate_col in enumerate(rate_cols, start=1):
                cov_col = companion_column(rate_col, "_coverage", fieldnames)
                alt_col = companion_column(rate_col, "_alt_count", fieldnames)
                record[f"T{i}_coverage"] = (
                    as_int(row.get(cov_col)) if cov_col else ""
                )
                record[f"T{i}_alt_count"] = (
                    as_int(row.get(alt_col)) if alt_col else ""
                )
                record[f"T{i}_edit_rate"] = rates[i - 1]

            record.update({
                "median_treated_edit_rate": statistics.median(rates),
                "mean_treated_edit_rate": statistics.mean(rates),
                "sd_treated_edit_rate": statistics.stdev(rates),
                "min_treated_edit_rate": min(rates),
                "max_treated_edit_rate": max(rates),
                "range_treated_edit_rate": max(rates) - min(rates),
                "label_recommendation": "median_treated_edit_rate",
                "label_caveat": (
                    "control_not_called; no background-corrected label "
                    "available from this table"
                ),
                "sequence_context": sequence,
                "sequence_orientation": orientation,
            })
            rows.append(record)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(output, "wt") as out:
        writer = csv.DictWriter(out, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    if fasta is not None:
        fasta.close()

    print(f"Rows written: {len(rows)}")
    print(f"Output: {output}")
    print(
        "Recommended provisional label: median_treated_edit_rate; "
        "do not treat missing control calls as zero."
    )


if __name__ == "__main__":
    main()
