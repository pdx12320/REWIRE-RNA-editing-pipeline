#!/usr/bin/env python3
"""Rebuild positive and strict-negative C sites from all covered transcript cytidines.

Every retained site has A/C/G/T depth strictly greater than the configured
threshold in every treated and control BAM. Reads are counted only when NH=1,
MAPQ and base quality pass, and duplicate/secondary/supplementary/QC-fail flags
are absent. The input universe is GTF-annotated exonic cytidines, independent of
whether a variant caller emitted the site.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import statistics
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


BASES = ("A", "C", "G", "T")
ATTRIBUTE_RE = re.compile(r'([^\s;]+)\s+"([^"]*)"')


def open_text(path: Path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def norm_chrom(chrom: str):
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def reverse_complement(sequence: str):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def parse_attributes(text: str):
    return dict(ATTRIBUTE_RE.findall(text))


@dataclass
class Exon:
    start0: int
    end0: int
    strand: str
    gene_id: str
    gene_name: str
    transcript_id: str


@dataclass
class Site:
    position: int
    ref: str
    alt: str
    strand: str
    gene_ids: set = field(default_factory=set)
    gene_names: set = field(default_factory=set)
    transcript_ids: set = field(default_factory=set)
    sample_counts: dict = field(default_factory=dict)


def read_manifest(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample", "group", "replicate"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Manifest must contain sample, group and replicate columns")
    if len({row["sample"] for row in rows}) != len(rows):
        raise ValueError("Manifest contains duplicate sample names")
    groups = defaultdict(list)
    for row in rows:
        groups[row["group"].lower()].append(row)
    if len(groups.get("treated", [])) != 3 or len(groups.get("control", [])) != 3:
        raise ValueError("Strict CU5.17 rebuild requires exactly 3 treated and 3 control samples")
    return rows


def read_gtf_exons(path: Path):
    exons = defaultdict(list)
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            chrom, _, _, start, end, _, strand, _, attributes_text = fields
            if strand not in {"+", "-"}:
                continue
            attributes = parse_attributes(attributes_text)
            gene_id = attributes.get("gene_id", "")
            if not gene_id:
                raise ValueError("GTF exon lacks gene_id at line {}".format(line_number))
            exons[chrom].append(
                Exon(
                    start0=int(start) - 1,
                    end0=int(end),
                    strand=strand,
                    gene_id=gene_id,
                    gene_name=attributes.get("gene_name", gene_id),
                    transcript_id=attributes.get("transcript_id", ""),
                )
            )
    if not exons:
        raise ValueError("GTF contains no exon features")
    return exons


def merge_intervals(exons):
    intervals = sorted((exon.start0, exon.end0) for exon in exons)
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def build_sites(reference, ref_chrom: str, exons, window_size: int = 101):
    sites = {}
    contig_length = reference.get_reference_length(ref_chrom)
    half = window_size // 2
    for exon in exons:
        start0 = max(exon.start0, 0)
        end0 = min(exon.end0, contig_length)
        if start0 >= end0:
            continue
        sequence = reference.fetch(ref_chrom, start0, end0).upper()
        wanted = "C" if exon.strand == "+" else "G"
        alt = "T" if exon.strand == "+" else "A"
        for offset, base in enumerate(sequence):
            if base != wanted:
                continue
            position = start0 + offset + 1
            position0 = position - 1
            # This strict implementation uses a contiguous genomic 101-mer for
            # BWA mappability. Keep only sites whose full window lies within an
            # annotated exon, so the window is also a valid transcript segment.
            if position0 - half < start0 or position0 + half >= end0:
                continue
            key = (position, wanted, alt)
            site = sites.get(key)
            if site is None:
                site = Site(position=position, ref=wanted, alt=alt, strand=exon.strand)
                sites[key] = site
            site.gene_ids.add(exon.gene_id)
            site.gene_names.add(exon.gene_name)
            if exon.transcript_id:
                site.transcript_ids.add(exon.transcript_id)
    return sites


def strict_read_counts(column, ref: str, alt: str, min_mapq: int, min_baseq: int):
    counts = {base: 0 for base in BASES}
    exclusions = defaultdict(int)
    for pileup_read in column.pileups:
        read = pileup_read.alignment
        if read.is_unmapped:
            continue
        if read.is_duplicate:
            exclusions["duplicate"] += 1
            continue
        if read.is_secondary:
            exclusions["secondary"] += 1
            continue
        if read.is_supplementary:
            exclusions["supplementary"] += 1
            continue
        if read.is_qcfail:
            exclusions["qcfail"] += 1
            continue
        if read.mapping_quality < min_mapq:
            exclusions["low_mapq"] += 1
            continue
        if not read.has_tag("NH"):
            exclusions["missing_nh"] += 1
            continue
        if int(read.get_tag("NH")) != 1:
            exclusions["nh_multimapper"] += 1
            continue
        if pileup_read.is_del or pileup_read.is_refskip or pileup_read.query_position is None:
            exclusions["deletion_or_refskip"] += 1
            continue
        query_position = pileup_read.query_position
        if read.query_qualities is None or read.query_qualities[query_position] < min_baseq:
            exclusions["low_baseq"] += 1
            continue
        if read.query_sequence is None:
            exclusions["non_acgt"] += 1
            continue
        base = read.query_sequence[query_position].upper()
        if base not in counts:
            exclusions["non_acgt"] += 1
            continue
        counts[base] += 1
    allele_depth = counts[ref] + counts[alt]
    return {
        "acgt_depth": sum(counts.values()),
        "ref_count": counts[ref],
        "alt_count": counts[alt],
        "allele_depth": allele_depth,
        "edit_rate": counts[alt] / allele_depth if allele_depth else None,
        "excluded_missing_nh": exclusions["missing_nh"],
        "excluded_nh_multimapper": exclusions["nh_multimapper"],
    }


def resolve_bam(bam_dir: Path, sample: str):
    path = bam_dir / "{}.splitncigarreads.bam".format(sample)
    if not path.is_file():
        raise FileNotFoundError("Missing uniform split BAM: {}".format(path))
    return path


def count_contig(alignment, bam_chrom: str, sites, intervals, min_mapq: int, min_baseq: int, max_depth: int):
    by_position = {site.position: site for site in sites.values()}
    observed = set()
    for start0, end0 in intervals:
        for column in alignment.pileup(
            bam_chrom,
            start0,
            end0,
            truncate=True,
            stepper="all",
            min_base_quality=0,
            min_mapping_quality=0,
            max_depth=max_depth,
        ):
            position = column.reference_pos + 1
            site = by_position.get(position)
            if site is None or position in observed:
                continue
            observed.add(position)
            site.sample_counts["__current__"] = strict_read_counts(
                column, site.ref, site.alt, min_mapq, min_baseq
            )
    for site in sites.values():
        site.sample_counts.setdefault(
            "__current__",
            {
                "acgt_depth": 0,
                "ref_count": 0,
                "alt_count": 0,
                "allele_depth": 0,
                "edit_rate": None,
                "excluded_missing_nh": 0,
                "excluded_nh_multimapper": 0,
            },
        )


def classify_site(site: Site, manifest, positive_threshold: float, min_positive_reps: int, min_alt_count: int, max_control_median: float):
    treated = [site.sample_counts[row["sample"]] for row in manifest if row["group"].lower() == "treated"]
    controls = [site.sample_counts[row["sample"]] for row in manifest if row["group"].lower() == "control"]
    treated_rates = [item["edit_rate"] for item in treated]
    control_rates = [item["edit_rate"] for item in controls]
    if any(value is None for value in treated_rates + control_rates):
        raise AssertionError("Fully covered site has a missing edit rate")
    treated_median = statistics.median(treated_rates)
    control_median = statistics.median(control_rates)
    raw = treated_median - control_median
    corrected = max(raw, 0.0)
    positive_reps = sum(item["alt_count"] >= min_alt_count for item in treated)
    strict_negative = all(item["alt_count"] == 0 for item in treated + controls)
    positive = (
        corrected >= positive_threshold
        and positive_reps >= min_positive_reps
        and control_median <= max_control_median
    )
    label_class = "positive" if positive else "strict_negative" if strict_negative else "intermediate"
    return {
        "treated_median": treated_median,
        "control_median": control_median,
        "raw_edit_rate_difference": raw,
        "corrected_editing_efficiency": corrected,
        "positive_treated_replicates": positive_reps,
        "label_class": label_class,
    }


def sequence_window(reference, chrom: str, site: Site, window_size: int):
    half = window_size // 2
    start0 = site.position - 1 - half
    end0 = site.position + half
    if start0 < 0 or end0 > reference.get_reference_length(chrom):
        return None
    sequence = reference.fetch(chrom, start0, end0).upper()
    if site.strand == "-":
        sequence = reverse_complement(sequence)
    if len(sequence) != window_size or sequence[half] != "C":
        raise AssertionError("Transcript-oriented window validation failed")
    return sequence


def format_number(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return "{:.10g}".format(value)
    return value


def sha256_file(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build(args):
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required") from exc

    manifest = read_manifest(args.manifest)
    if args.expected_samples and len(manifest) != args.expected_samples:
        raise ValueError("Expected {} samples, found {}".format(args.expected_samples, len(manifest)))
    exons_by_chrom = read_gtf_exons(args.gtf)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "status": "running",
        "min_acgt_depth_exclusive": args.min_depth,
        "positive_threshold": args.positive_threshold,
        "strict_negative_rule": "alt_count == 0 in all six samples",
        "candidate_universe": "all GTF-exonic transcript cytidines with a contiguous 101-nt exon window",
        "read_rule": "primary, nonduplicate, non-QC-fail, MAPQ>=threshold, baseQ>=threshold, NH=1",
        "samples": [row["sample"] for row in manifest],
        "contigs": {},
        "rows": 0,
        "positive": 0,
        "strict_negative": 0,
        "intermediate": 0,
    }
    sample_fields = []
    for row in manifest:
        sample = row["sample"]
        sample_fields.extend(
            [
                sample + "_acgt_depth",
                sample + "_ref_count",
                sample + "_alt_count",
                sample + "_allele_depth",
                sample + "_edit_rate",
                sample + "_excluded_missing_nh",
                sample + "_excluded_nh_multimapper",
            ]
        )
    fields = [
        "chrom", "position", "ref", "alt", "strand", "gene_id", "gene_name", "transcript_id",
        "gene_ambiguous", "sequence_context", "sequence_length", "center_index", "center_base",
        "orientation_qc", "transcript_oriented_ref", "transcript_oriented_alt",
        "minimum_acgt_depth", "sufficient_coverage_all_six",
    ] + sample_fields + [
        "treated_median", "control_median", "raw_edit_rate_difference",
        "corrected_editing_efficiency", "positive_treated_replicates", "label_class",
        "training_eligible", "label_confidence", "elevated_control_background", "exclusion_reason",
    ]

    with pysam.FastaFile(str(args.reference)) as reference, open_text(args.output, "wt") as output:
        ref_map = {norm_chrom(name): name for name in reference.references}
        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gtf_chrom in sorted(exons_by_chrom, key=lambda value: (norm_chrom(value), value)):
            ref_chrom = ref_map.get(norm_chrom(gtf_chrom))
            if ref_chrom is None:
                continue
            exons = exons_by_chrom[gtf_chrom]
            sites = build_sites(reference, ref_chrom, exons, args.window_size)
            initial_count = len(sites)
            intervals = merge_intervals(exons)
            for meta in manifest:
                bam = resolve_bam(args.bam_dir, meta["sample"])
                with pysam.AlignmentFile(str(bam), "rb") as alignment:
                    bam_map = {norm_chrom(name): name for name in alignment.references}
                    bam_chrom = bam_map.get(norm_chrom(gtf_chrom))
                    if bam_chrom is None:
                        raise ValueError("Contig {} is absent from {}".format(gtf_chrom, bam))
                    count_contig(
                        alignment,
                        bam_chrom,
                        sites,
                        intervals,
                        args.min_mapq,
                        args.min_baseq,
                        args.max_depth,
                    )
                retained = {}
                for key, site in sites.items():
                    counts = site.sample_counts.pop("__current__")
                    if counts["acgt_depth"] > args.min_depth:
                        site.sample_counts[meta["sample"]] = counts
                        retained[key] = site
                sites = retained
                if not sites:
                    break

            emitted = 0
            for site in sites.values():
                sequence = sequence_window(reference, ref_chrom, site, args.window_size)
                if sequence is None:
                    continue
                labels = classify_site(
                    site,
                    manifest,
                    args.positive_threshold,
                    args.min_positive_treated_reps,
                    args.min_positive_alt_count,
                    args.max_positive_control_median,
                )
                gene_ids = sorted(site.gene_ids)
                gene_names = sorted(site.gene_names)
                gene_ambiguous = len(gene_ids) != 1
                row = {
                    "chrom": ref_chrom,
                    "position": site.position,
                    "ref": site.ref,
                    "alt": site.alt,
                    "strand": site.strand,
                    "gene_id": "|".join(gene_ids),
                    "gene_name": "|".join(gene_names),
                    "transcript_id": "|".join(sorted(site.transcript_ids)),
                    "gene_ambiguous": int(gene_ambiguous),
                    "sequence_context": sequence,
                    "sequence_length": args.window_size,
                    "center_index": args.window_size // 2,
                    "center_base": "C",
                    "orientation_qc": "pass",
                    "transcript_oriented_ref": "C",
                    "transcript_oriented_alt": "T",
                    "minimum_acgt_depth": min(
                        site.sample_counts[meta["sample"]]["acgt_depth"] for meta in manifest
                    ),
                    "sufficient_coverage_all_six": 1,
                    **labels,
                    "training_eligible": int(not gene_ambiguous and labels["label_class"] != "intermediate"),
                    "label_confidence": "high",
                    "elevated_control_background": int(labels["control_median"] > args.max_positive_control_median),
                    "exclusion_reason": (
                        "ambiguous_gene" if gene_ambiguous else "intermediate_editing" if labels["label_class"] == "intermediate" else "none"
                    ),
                }
                for meta in manifest:
                    sample = meta["sample"]
                    counts = site.sample_counts[sample]
                    for field_name, value in counts.items():
                        row[sample + "_" + field_name] = format_number(value)
                writer.writerow({field: format_number(row.get(field, "NA")) for field in fields})
                emitted += 1
                summary["rows"] += 1
                summary[labels["label_class"]] += 1
            summary["contigs"][ref_chrom] = {
                "annotated_cytidines_with_contiguous_101nt_exon_window": initial_count,
                "emitted_full_coverage": emitted,
            }

    summary["status"] = "pass"
    summary["output"] = str(args.output.resolve())
    summary["output_sha256"] = sha256_file(args.output)
    audit_path = Path(str(args.output) + ".audit.json")
    audit_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--bam-dir", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-samples", type=int, default=6)
    parser.add_argument("--min-depth", type=int, default=50, help="Exclusive threshold; default retains depth >50")
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--max-depth", type=int, default=1000000)
    parser.add_argument("--window-size", type=int, default=101)
    parser.add_argument("--positive-threshold", type=float, default=0.10)
    parser.add_argument("--min-positive-treated-reps", type=int, default=2)
    parser.add_argument("--min-positive-alt-count", type=int, default=3)
    parser.add_argument("--max-positive-control-median", type=float, default=0.01)
    args = parser.parse_args(argv)
    if args.window_size != 101:
        parser.error("This LAMAR handoff is fixed to 101 nt")
    if args.min_depth < 0 or args.max_depth < 1:
        parser.error("Depth arguments are invalid")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        summary = build(args)
    except Exception as error:
        raise SystemExit("ERROR: {}".format(error))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
