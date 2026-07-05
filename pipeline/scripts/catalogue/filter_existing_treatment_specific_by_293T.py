#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def norm_chrom(chrom):
    chrom = chrom.strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def find_column(fieldnames, candidates, label):
    lookup = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    raise SystemExit(
        f"ERROR: cannot find {label} column. "
        f"Available columns: {', '.join(fieldnames)}"
    )


def load_catalogue(path):
    variants = set()
    with open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            filt = fields[6] if len(fields) > 6 else "."
            if filt not in {".", "PASS"}:
                continue
            chrom = norm_chrom(fields[0])
            pos = int(fields[1])
            ref = fields[3].upper()
            for alt in fields[4].split(","):
                variants.add((chrom, pos, ref, alt.upper()))
    return variants


def write_rows(path, fieldnames, rows):
    with open_text(path, "wt") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter an already completed REWIRE treatment-specific table "
            "against the GRCh38 293T genomic catalogue by exact "
            "CHROM:POS:REF:ALT matching."
        )
    )
    parser.add_argument("--treatment-specific", required=True)
    parser.add_argument("--catalogue-vcf", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    input_table = Path(args.treatment_specific)
    catalogue_vcf = Path(args.catalogue_vcf)
    output_dir = Path(args.output_dir)

    if not input_table.is_file():
        raise SystemExit(f"ERROR: input table not found: {input_table}")
    if not catalogue_vcf.is_file():
        raise SystemExit(f"ERROR: catalogue VCF not found: {catalogue_vcf}")

    print(f"Loading catalogue: {catalogue_vcf}")
    catalogue = load_catalogue(catalogue_vcf)
    print(f"Catalogue exact alleles: {len(catalogue)}")

    with open_text(input_table) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: treatment-specific table has no header")

        original_fields = list(reader.fieldnames)
        chrom_col = find_column(
            original_fields,
            ["chrom", "chr", "chromosome", "region"],
            "chromosome",
        )
        pos_col = find_column(
            original_fields,
            ["position", "pos", "coordinate"],
            "position",
        )
        ref_col = find_column(
            original_fields,
            ["ref", "reference"],
            "reference allele",
        )
        alt_col = find_column(
            original_fields,
            ["alt", "alternate", "alternative"],
            "alternate allele",
        )

        output_fields = list(original_fields)
        if "genomic_catalogue_overlap" not in output_fields:
            output_fields.append("genomic_catalogue_overlap")

        all_rows = []
        retained = []
        excluded = []

        for row in reader:
            key = (
                norm_chrom(row[chrom_col]),
                int(row[pos_col]),
                row[ref_col].upper(),
                row[alt_col].upper(),
            )
            overlap = int(key in catalogue)
            row["genomic_catalogue_overlap"] = overlap
            all_rows.append(row)
            if overlap:
                excluded.append(row)
            else:
                retained.append(row)

    output_dir.mkdir(parents=True, exist_ok=True)

    annotated_path = output_dir / "CU5.17_EGFP_GC.treatment_specific.before_catalogue.annotated.tsv.gz"
    final_path = output_dir / "CU5.17_EGFP_GC.treatment_specific.tsv.gz"
    excluded_path = output_dir / "CU5.17_EGFP_GC.catalogue_overlaps.tsv.gz"
    summary_path = output_dir / "catalogue_integration_summary.tsv"

    write_rows(annotated_path, output_fields, all_rows)
    write_rows(final_path, output_fields, retained)
    write_rows(excluded_path, output_fields, excluded)

    with open(summary_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["metric", "count"])
        writer.writerow(["treatment_specific_before_catalogue", len(all_rows)])
        writer.writerow(["catalogue_overlap", len(excluded)])
        writer.writerow(["final_treatment_specific", len(retained)])
        writer.writerow(["catalogue_exact_alleles", len(catalogue)])

    print(f"treatment_specific_before_catalogue={len(all_rows)}")
    print(f"catalogue_overlap={len(excluded)}")
    print(f"final_treatment_specific={len(retained)}")
    print(f"Final table: {final_path}")
    print(f"Excluded overlaps: {excluded_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
