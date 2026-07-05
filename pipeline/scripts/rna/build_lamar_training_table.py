#!/usr/bin/env python3
"""Build a replicate-aware LAMAR training table from candidate-site base counts.

The output keeps all per-sample counts, creates continuous treated/control editing
labels, reports a pooled Fisher screening statistic, and emits a transcript-oriented
sequence window centered on the editable cytidine.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import statistics
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


def open_text(path: Path | str, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def norm_chrom(chrom: str) -> str:
    chrom = chrom.strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def canonical_key(chrom: str, position: str | int, ref: str, alt: str) -> Tuple[str, int, str, str]:
    return norm_chrom(chrom), int(position), ref.upper(), alt.upper()


def find_column(fieldnames: Sequence[str], candidates: Sequence[str], label: str, required: bool = True):
    lookup = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    if required:
        raise SystemExit(
            f"ERROR: cannot find {label} column. Available columns: {', '.join(fieldnames)}"
        )
    return None


def read_manifest(path: Path) -> List[dict]:
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample", "group", "replicate"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("ERROR: manifest must contain sample, group and replicate columns")
    groups = {row["group"] for row in rows}
    if not {"treated", "control"}.issubset(groups):
        raise SystemExit("ERROR: manifest must include treated and control samples")
    return rows


def read_keyed_table(path: Path) -> Tuple[List[str], OrderedDict]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: table has no header: {path}")
        fields = list(reader.fieldnames)
        chrom_col = find_column(fields, ["chrom", "chr", "chromosome", "region"], "chromosome")
        pos_col = find_column(fields, ["position", "pos", "coordinate"], "position")
        ref_col = find_column(fields, ["ref", "reference"], "reference allele")
        alt_col = find_column(fields, ["alt", "alternate", "alternative"], "alternate allele")
        rows = OrderedDict()
        for row in reader:
            key = canonical_key(row[chrom_col], row[pos_col], row[ref_col], row[alt_col])
            if key in rows:
                previous = rows[key]
                previous_strand = previous.get("vep_strand", previous.get("strand", ""))
                current_strand = row.get("vep_strand", row.get("strand", ""))
                if previous_strand and current_strand and previous_strand != current_strand:
                    raise SystemExit(f"ERROR: conflicting strand annotation for {key}")
                continue
            rows[key] = row
    return fields, rows


def read_count_table(path: Path) -> Dict[Tuple[str, int, str, str], dict]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: count table has no header: {path}")
        required = {"chrom", "position", "ref", "alt", "ref_count", "alt_count", "allele_depth", "acgt_depth"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise SystemExit(f"ERROR: {path} missing columns: {', '.join(sorted(missing))}")
        return {
            canonical_key(row["chrom"], row["position"], row["ref"], row["alt"]): row
            for row in reader
        }


def read_metadata(path: Path | None):
    if path is None:
        return [], {}
    fields, rows = read_keyed_table(path)
    key_aliases = {
        "chrom", "chr", "chromosome", "region", "position", "pos", "coordinate",
        "ref", "reference", "alt", "alternate", "alternative",
    }
    extra_fields = [field for field in fields if field.lower() not in key_aliases]
    return extra_fields, rows


def parse_int(value, default=0) -> int:
    if value in (None, "", "NA", "."):
        return default
    return int(float(value))


def median(values: Sequence[float]) -> float:
    return statistics.median(values) if values else float("nan")


def median_abs_deviation(values: Sequence[float]) -> float:
    if not values:
        return float("nan")
    center = statistics.median(values)
    return statistics.median(abs(value - center) for value in values)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for [[a, b], [c, d]]."""
    if min(a, b, c, d) < 0:
        raise ValueError("Fisher table counts must be non-negative")
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2
    if total == 0:
        return 1.0

    def log_choose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return float("-inf")
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def probability(x: int) -> float:
        return math.exp(log_choose(col1, x) + log_choose(total - col1, row1 - x) - log_choose(total, row1))

    lower = max(0, row1 - (total - col1))
    upper = min(row1, col1)
    observed = probability(a)
    p_value = 0.0
    for x in range(lower, upper + 1):
        candidate = probability(x)
        if candidate <= observed * (1 + 1e-12):
            p_value += candidate
    return min(1.0, p_value)


def bh_fdr(p_values: Sequence[float]) -> List[float]:
    if not p_values:
        return []
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_from_end, (index, p_value) in enumerate(reversed(indexed), start=1):
        rank = total - rank_from_end + 1
        running = min(running, p_value * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def fasta_name_map(references: Sequence[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for reference in references:
        key = norm_chrom(reference)
        if key in mapping and mapping[key] != reference:
            raise SystemExit(
                f"ERROR: ambiguous FASTA contigs after chr normalization: {mapping[key]} and {reference}"
            )
        mapping[key] = reference
    return mapping


def fetch_oriented_window(fasta, chrom_map: Mapping[str, str], key, strand: int, window_size: int) -> str:
    chrom, position, _ref, _alt = key
    fasta_chrom = chrom_map.get(norm_chrom(chrom))
    if fasta_chrom is None:
        return "N" * window_size
    half = window_size // 2
    center0 = position - 1
    start = center0 - half
    end = center0 + half + 1
    reference_length = fasta.get_reference_length(fasta_chrom)
    left_pad = max(0, -start)
    right_pad = max(0, end - reference_length)
    sequence = (
        "N" * left_pad
        + fasta.fetch(fasta_chrom, max(0, start), min(end, reference_length)).upper()
        + "N" * right_pad
    )
    if len(sequence) != window_size:
        sequence = (sequence + "N" * window_size)[:window_size]
    return reverse_complement(sequence) if strand == -1 else sequence


def format_number(value, digits: int = 10):
    if isinstance(value, float) and math.isnan(value):
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return value


def infer_strand(row: Mapping[str, str], key) -> int:
    strand_value = row.get("vep_strand", row.get("strand", ""))
    try:
        strand = int(float(strand_value))
    except (TypeError, ValueError):
        _chrom, _position, ref, alt = key
        if ref == "C" and alt == "T":
            strand = 1
        elif ref == "G" and alt == "A":
            strand = -1
        else:
            strand = 0
    return strand if strand in (-1, 1) else 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine six-sample base counts into continuous, background-corrected LAMAR labels."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--candidate-table", required=True, type=Path)
    parser.add_argument("--count-dir", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--site-metadata", type=Path, default=None, help="Optional exact-allele TSV with PUF/site features")
    parser.add_argument("--window-size", type=int, default=101)
    parser.add_argument("--min-allele-depth", type=int, default=20)
    parser.add_argument("--positive-threshold", type=float, default=0.01)
    parser.add_argument("--near-background-threshold", type=float, default=0.002)
    args = parser.parse_args()

    if args.window_size < 3 or args.window_size % 2 != 1:
        parser.error("window-size must be an odd integer of at least 3")
    if args.min_allele_depth < 1:
        parser.error("min-allele-depth must be positive")
    if not 0 <= args.near_background_threshold <= args.positive_threshold <= 1:
        parser.error("require 0 <= near-background-threshold <= positive-threshold <= 1")

    try:
        import pysam
    except ImportError as exc:
        raise SystemExit(
            "ERROR: pysam is required. Create the environment with pipeline/env/lamar_labels.yml"
        ) from exc

    manifest = read_manifest(args.manifest)
    original_fields, candidates = read_keyed_table(args.candidate_table)
    metadata_fields, metadata = read_metadata(args.site_metadata)
    sample_counts = {}
    for meta in manifest:
        path = args.count_dir / f"{meta['sample']}.candidate_base_counts.tsv.gz"
        if not path.is_file():
            raise SystemExit(f"ERROR: missing count table: {path}")
        sample_counts[meta["sample"]] = read_count_table(path)

    for sample, table in sample_counts.items():
        missing = [key for key in candidates if key not in table]
        if missing:
            raise SystemExit(f"ERROR: {sample} count table is missing {len(missing)} candidate sites")

    treated_samples = [row["sample"] for row in manifest if row["group"] == "treated"]
    control_samples = [row["sample"] for row in manifest if row["group"] == "control"]

    derived_fields = [
        "transcript_oriented_sequence", "sequence_window_size", "editable_base_index_0based",
        "sequence_center_base", "sequence_center_is_c", "strand_rule_valid",
        "treated_covered_reps", "control_covered_reps", "all_replicates_depth_pass",
        "minimum_allele_depth", "treated_median_edit_rate", "control_median_edit_rate",
        "raw_edit_rate_difference", "corrected_editing_efficiency", "treated_pooled_edit_rate",
        "control_pooled_edit_rate", "treated_edit_rate_mad", "control_edit_rate_mad",
        "pooled_log2_odds_ratio", "pooled_fisher_p", "pooled_fisher_fdr",
        "label_class", "label_confidence", "training_eligible",
    ]
    sample_fields: List[str] = []
    for meta in manifest:
        sample = meta["sample"]
        sample_fields += [
            f"{sample}_acgt_depth", f"{sample}_ref_count", f"{sample}_alt_count",
            f"{sample}_allele_depth", f"{sample}_edit_rate", f"{sample}_alt_forward_count",
            f"{sample}_alt_reverse_count", f"{sample}_alt_strand_balance",
        ]

    output_fields = list(original_fields)
    for field in metadata_fields + sample_fields + derived_fields:
        if field not in output_fields:
            output_fields.append(field)

    rows: List[dict] = []
    p_values: List[float] = []
    with pysam.FastaFile(str(args.reference)) as fasta:
        chrom_map = fasta_name_map(fasta.references)
        for key, original in candidates.items():
            row = dict(original)
            if key in metadata:
                for field in metadata_fields:
                    row[field] = metadata[key].get(field, "")

            strand = infer_strand(row, key)
            _chrom, _position, ref, alt = key
            strand_rule_valid = int(
                (strand == 1 and ref == "C" and alt == "T")
                or (strand == -1 and ref == "G" and alt == "A")
            )
            sequence = fetch_oriented_window(fasta, chrom_map, key, strand, args.window_size)
            center = sequence[args.window_size // 2]

            group_rates = {"treated": [], "control": []}
            group_ref = {"treated": 0, "control": 0}
            group_alt = {"treated": 0, "control": 0}
            group_covered = {"treated": 0, "control": 0}
            minimum_depth = None

            for meta in manifest:
                sample = meta["sample"]
                group = meta["group"]
                count = sample_counts[sample][key]
                ref_count = parse_int(count["ref_count"])
                alt_count = parse_int(count["alt_count"])
                allele_depth = parse_int(count["allele_depth"])
                acgt_depth = parse_int(count["acgt_depth"])
                rate = alt_count / allele_depth if allele_depth else float("nan")
                minimum_depth = allele_depth if minimum_depth is None else min(minimum_depth, allele_depth)
                if allele_depth >= args.min_allele_depth:
                    group_covered[group] += 1
                    group_rates[group].append(rate)
                    group_ref[group] += ref_count
                    group_alt[group] += alt_count

                row[f"{sample}_acgt_depth"] = acgt_depth
                row[f"{sample}_ref_count"] = ref_count
                row[f"{sample}_alt_count"] = alt_count
                row[f"{sample}_allele_depth"] = allele_depth
                row[f"{sample}_edit_rate"] = format_number(rate)
                row[f"{sample}_alt_forward_count"] = parse_int(count.get("alt_forward_count"))
                row[f"{sample}_alt_reverse_count"] = parse_int(count.get("alt_reverse_count"))
                row[f"{sample}_alt_strand_balance"] = count.get("alt_strand_balance", "NA")

            treated_median = median(group_rates["treated"])
            control_median = median(group_rates["control"])
            raw_difference = (
                treated_median - control_median
                if not math.isnan(treated_median) and not math.isnan(control_median)
                else float("nan")
            )
            corrected = max(0.0, raw_difference) if not math.isnan(raw_difference) else float("nan")
            treated_denominator = group_ref["treated"] + group_alt["treated"]
            control_denominator = group_ref["control"] + group_alt["control"]
            treated_pooled = group_alt["treated"] / treated_denominator if treated_denominator else float("nan")
            control_pooled = group_alt["control"] / control_denominator if control_denominator else float("nan")
            fisher_p = fisher_exact_two_sided(
                group_alt["treated"], group_ref["treated"], group_alt["control"], group_ref["control"]
            )
            p_values.append(fisher_p)
            odds_ratio = (
                (group_alt["treated"] + 0.5) * (group_ref["control"] + 0.5)
                / ((group_ref["treated"] + 0.5) * (group_alt["control"] + 0.5))
            )

            all_depth_pass = int(
                group_covered["treated"] == len(treated_samples)
                and group_covered["control"] == len(control_samples)
            )
            catalogue_overlap = parse_int(row.get("genomic_catalogue_overlap"), 0)
            sequence_center_is_c = int(center == "C")
            training_eligible = int(
                all_depth_pass == 1 and strand_rule_valid == 1 and sequence_center_is_c == 1 and catalogue_overlap == 0
            )

            if math.isnan(corrected):
                label_class = "unresolved"
            elif corrected >= args.positive_threshold:
                label_class = "positive"
            elif corrected <= args.near_background_threshold:
                label_class = "near_background"
            else:
                label_class = "intermediate"

            if training_eligible:
                label_confidence = "high"
            elif (
                group_covered["treated"] >= max(1, len(treated_samples) - 1)
                and group_covered["control"] >= max(1, len(control_samples) - 1)
                and strand_rule_valid and sequence_center_is_c and not catalogue_overlap
            ):
                label_confidence = "moderate"
            else:
                label_confidence = "low"

            row.update(
                {
                    "transcript_oriented_sequence": sequence,
                    "sequence_window_size": args.window_size,
                    "editable_base_index_0based": args.window_size // 2,
                    "sequence_center_base": center,
                    "sequence_center_is_c": sequence_center_is_c,
                    "strand_rule_valid": strand_rule_valid,
                    "treated_covered_reps": group_covered["treated"],
                    "control_covered_reps": group_covered["control"],
                    "all_replicates_depth_pass": all_depth_pass,
                    "minimum_allele_depth": minimum_depth if minimum_depth is not None else 0,
                    "treated_median_edit_rate": format_number(treated_median),
                    "control_median_edit_rate": format_number(control_median),
                    "raw_edit_rate_difference": format_number(raw_difference),
                    "corrected_editing_efficiency": format_number(corrected),
                    "treated_pooled_edit_rate": format_number(treated_pooled),
                    "control_pooled_edit_rate": format_number(control_pooled),
                    "treated_edit_rate_mad": format_number(median_abs_deviation(group_rates["treated"])),
                    "control_edit_rate_mad": format_number(median_abs_deviation(group_rates["control"])),
                    "pooled_log2_odds_ratio": format_number(math.log2(odds_ratio)),
                    "pooled_fisher_p": format_number(fisher_p),
                    "pooled_fisher_fdr": "NA",
                    "label_class": label_class,
                    "label_confidence": label_confidence,
                    "training_eligible": training_eligible,
                }
            )
            rows.append(row)

    for row, fdr in zip(rows, bh_fdr(p_values)):
        row["pooled_fisher_fdr"] = format_number(fdr)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output, "wt") as output:
        writer = csv.DictWriter(output, fieldnames=output_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    eligible = sum(parse_int(row["training_eligible"]) for row in rows)
    high = sum(row["label_confidence"] == "high" for row in rows)
    positive = sum(row["label_class"] == "positive" for row in rows)
    print(f"training rows={len(rows)}")
    print(f"training eligible={eligible}")
    print(f"high-confidence labels={high}")
    print(f"positive labels={positive}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
