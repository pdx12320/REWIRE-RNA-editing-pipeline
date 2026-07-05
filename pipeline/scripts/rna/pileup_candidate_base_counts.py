#!/usr/bin/env python3
"""Measure quality-filtered A/C/G/T counts at candidate sites in every sample.

Unlike a caller table, this script records counts even when the alternate allele does
not cross a calling threshold. This is required for continuous editing-efficiency
labels and for distinguishing a true zero from a control non-call.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

BASES = ("A", "C", "G", "T")


def open_text(path: Path | str, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def norm_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def find_column(fieldnames: Sequence[str], candidates: Sequence[str], label: str) -> str:
    lookup = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    raise SystemExit(
        f"ERROR: cannot find {label} column. Available columns: {', '.join(fieldnames)}"
    )


@dataclass(frozen=True, order=True)
class Site:
    chrom: str
    position: int
    ref: str
    alt: str


def read_manifest(path: Path) -> List[dict]:
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"sample", "group", "replicate"}.issubset(rows[0]):
        raise SystemExit("ERROR: manifest must contain sample, group and replicate columns")
    return rows


def read_sites(path: Path) -> List[Site]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: candidate table has no header: {path}")
        chrom_col = find_column(reader.fieldnames, ["chrom", "chr", "chromosome", "region"], "chromosome")
        pos_col = find_column(reader.fieldnames, ["position", "pos", "coordinate"], "position")
        ref_col = find_column(reader.fieldnames, ["ref", "reference"], "reference allele")
        alt_col = find_column(reader.fieldnames, ["alt", "alternate", "alternative"], "alternate allele")
        sites = {
            Site(
                chrom=row[chrom_col].strip(),
                position=int(row[pos_col]),
                ref=row[ref_col].upper(),
                alt=row[alt_col].upper(),
            )
            for row in reader
            if row.get(chrom_col) and row.get(pos_col) and row.get(ref_col) and row.get(alt_col)
        }
    return sorted(sites, key=lambda s: (norm_chrom(s.chrom), s.position, s.ref, s.alt))


def resolve_bam_path(bam_dir: Path, sample: str, srr: str) -> Path:
    candidates = [
        bam_dir / f"{sample}.splitncigarreads.bam",
        bam_dir / f"{srr}.splitncigarreads.bam" if srr else Path("/__missing__"),
        bam_dir / f"{sample}.bam",
        bam_dir / f"{srr}.bam" if srr else Path("/__missing__"),
    ]
    for path in candidates:
        if path.is_file():
            return path
    checked = ", ".join(str(p) for p in candidates if str(p) != "/__missing__")
    raise SystemExit(f"ERROR: BAM not found for {sample}. Checked: {checked}")


def reference_name_map(references: Sequence[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for reference in references:
        key = norm_chrom(reference)
        if key in mapping and mapping[key] != reference:
            raise SystemExit(
                f"ERROR: ambiguous BAM contig names after chr normalization: {mapping[key]} and {reference}"
            )
        mapping[key] = reference
    return mapping


def count_site(alignment, bam_chrom: str, position: int, min_mapq: int, min_baseq: int, max_depth: int) -> dict:
    counts = {base: 0 for base in BASES}
    forward = {base: 0 for base in BASES}
    reverse = {base: 0 for base in BASES}
    excluded = {
        "low_mapq": 0,
        "low_baseq": 0,
        "duplicate": 0,
        "secondary": 0,
        "supplementary": 0,
        "qcfail": 0,
        "deletion_or_refskip": 0,
        "non_acgt": 0,
    }

    for column in alignment.pileup(
        bam_chrom,
        position - 1,
        position,
        truncate=True,
        stepper="all",
        min_base_quality=0,
        min_mapping_quality=0,
        max_depth=max_depth,
    ):
        if column.reference_pos != position - 1:
            continue
        for pileup_read in column.pileups:
            read = pileup_read.alignment
            if read.is_unmapped:
                continue
            if read.is_duplicate:
                excluded["duplicate"] += 1
                continue
            if read.is_secondary:
                excluded["secondary"] += 1
                continue
            if read.is_supplementary:
                excluded["supplementary"] += 1
                continue
            if read.is_qcfail:
                excluded["qcfail"] += 1
                continue
            if read.mapping_quality < min_mapq:
                excluded["low_mapq"] += 1
                continue
            if pileup_read.is_del or pileup_read.is_refskip or pileup_read.query_position is None:
                excluded["deletion_or_refskip"] += 1
                continue
            query_position = pileup_read.query_position
            qualities = read.query_qualities
            if qualities is None or qualities[query_position] < min_baseq:
                excluded["low_baseq"] += 1
                continue
            sequence = read.query_sequence
            if sequence is None:
                excluded["non_acgt"] += 1
                continue
            base = sequence[query_position].upper()
            if base not in counts:
                excluded["non_acgt"] += 1
                continue
            counts[base] += 1
            (reverse if read.is_reverse else forward)[base] += 1

    result = {}
    for base in BASES:
        result[f"{base}_count"] = counts[base]
        result[f"{base}_forward_count"] = forward[base]
        result[f"{base}_reverse_count"] = reverse[base]
    result["acgt_depth"] = sum(counts.values())
    result["excluded_read_observations"] = sum(excluded.values())
    for reason, value in excluded.items():
        result[f"excluded_{reason}"] = value
    return result


def output_fields() -> List[str]:
    fields = ["sample", "group", "replicate", "chrom", "position", "ref", "alt"]
    for base in BASES:
        fields += [f"{base}_count", f"{base}_forward_count", f"{base}_reverse_count"]
    fields += [
        "acgt_depth",
        "ref_count",
        "alt_count",
        "allele_depth",
        "edit_rate",
        "alt_forward_count",
        "alt_reverse_count",
        "alt_strand_balance",
        "excluded_read_observations",
        "excluded_low_mapq",
        "excluded_low_baseq",
        "excluded_duplicate",
        "excluded_secondary",
        "excluded_supplementary",
        "excluded_qcfail",
        "excluded_deletion_or_refskip",
        "excluded_non_acgt",
    ]
    return fields


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract quality-filtered A/C/G/T counts at every candidate site from all BAMs. "
            "Counts are emitted even for sites absent from a variant-caller output."
        )
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--candidate-table", required=True, type=Path)
    parser.add_argument("--bam-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--reference", type=Path, default=None, help="Reference FASTA; required for CRAM input")
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--max-depth", type=int, default=1_000_000)
    args = parser.parse_args()

    if args.min_mapq < 0 or args.min_baseq < 0 or args.max_depth < 1:
        parser.error("quality thresholds must be non-negative and max-depth must be positive")

    try:
        import pysam
    except ImportError as exc:
        raise SystemExit(
            "ERROR: pysam is required. Create the environment with pipeline/env/lamar_labels.yml"
        ) from exc

    manifest = read_manifest(args.manifest)
    sites = read_sites(args.candidate_table)
    if not sites:
        raise SystemExit("ERROR: candidate table contains no sites")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = output_fields()
    reference_filename = str(args.reference) if args.reference else None

    for meta in manifest:
        sample = meta["sample"]
        bam_path = resolve_bam_path(args.bam_dir, sample, meta.get("srr", ""))
        output_path = args.output_dir / f"{sample}.candidate_base_counts.tsv.gz"
        mode = "rc" if bam_path.suffix.lower() == ".cram" else "rb"
        with pysam.AlignmentFile(bam_path, mode, reference_filename=reference_filename) as alignment, open_text(
            output_path, "wt"
        ) as output:
            chrom_map = reference_name_map(alignment.references)
            writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for site in sites:
                bam_chrom = chrom_map.get(norm_chrom(site.chrom))
                if bam_chrom is None:
                    raise SystemExit(f"ERROR: candidate contig {site.chrom!r} is absent from BAM {bam_path}")
                measured = count_site(alignment, bam_chrom, site.position, args.min_mapq, args.min_baseq, args.max_depth)
                ref_count = measured.get(f"{site.ref}_count", 0)
                alt_count = measured.get(f"{site.alt}_count", 0)
                allele_depth = ref_count + alt_count
                alt_forward = measured.get(f"{site.alt}_forward_count", 0)
                alt_reverse = measured.get(f"{site.alt}_reverse_count", 0)
                row = {
                    "sample": sample,
                    "group": meta["group"],
                    "replicate": meta["replicate"],
                    "chrom": site.chrom,
                    "position": site.position,
                    "ref": site.ref,
                    "alt": site.alt,
                    **measured,
                    "ref_count": ref_count,
                    "alt_count": alt_count,
                    "allele_depth": allele_depth,
                    "edit_rate": f"{alt_count / allele_depth:.10g}" if allele_depth else "NA",
                    "alt_forward_count": alt_forward,
                    "alt_reverse_count": alt_reverse,
                    "alt_strand_balance": (
                        f"{min(alt_forward, alt_reverse) / max(alt_forward, alt_reverse):.10g}"
                        if max(alt_forward, alt_reverse) > 0
                        else "NA"
                    ),
                }
                writer.writerow(row)
        print(f"candidate base counts: {sample}\t{len(sites)}\t{output_path}")


if __name__ == "__main__":
    main()
