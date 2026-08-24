#!/usr/bin/env python3
"""Build a leakage-resistant, auditable LAMAR scalar-regression handoff.

This script consumes the frozen broad-matrix label and sequence-context tables.
It does not read BAM files or recompute scientific labels.  The split table
contains every training-eligible row; downstream tools select the recommended
high-confidence subset without changing the shared split assignment.
"""

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import random
import shutil
import sys
from collections import defaultdict
from pathlib import Path


KEY_FIELDS = ("chrom", "position", "ref", "alt")
SEQUENCE_LENGTH = 101
CENTER_INDEX = 50
DEFAULT_SEED = 20260715
SPLIT_FRACTIONS = {"train": 0.8, "validation": 0.1, "test": 0.1}
MISSING = {"", ".", "NA", "N/A", "nan", "NaN", "None", "null"}
EXPECTED_OUTPUTS = (
    "CU5.17_lamar_all_eligible.tsv.gz",
    "CU5.17_lamar_high_confidence.tsv.gz",
    "CU5.17_lamar_high_confidence_low_control.tsv.gz",
    "CU5.17_lamar_excluded.tsv.gz",
    "CU5.17_lamar_splits.tsv.gz",
    "data_dictionary.tsv",
    "split_qc.json",
    "handoff_manifest.json",
    "README.md",
    "checksums.sha256",
)
PUBLIC_OUTPUTS = (
    "CU5.17_lamar_all_eligible.tsv.gz",
    "CU5.17_lamar_high_confidence.tsv.gz",
    "CU5.17_lamar_splits.tsv.gz",
    "data_dictionary.tsv",
    "handoff_manifest.json",
)


def open_text(path, mode="rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_tsv(path):
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("TSV has no header: {}".format(path))
        rows = [dict(row) for row in reader]
    return list(reader.fieldnames), rows


def write_tsv(path, fields, rows):
    """Write stable TSV or gzip-compressed TSV (gzip mtime is fixed at zero)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = None
    if path.suffix == ".gz":
        raw = path.open("wb")
        compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        handle = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    else:
        handle = path.open("w", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_output(row.get(field)) for field in fields})
    finally:
        handle.close()
        if raw is not None:
            raw.close()


def format_output(value):
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "1" if value else "0"
    return str(value)


def allele_key(row):
    try:
        position = int(row["position"])
    except (KeyError, TypeError, ValueError):
        raise ValueError("Invalid 1-based position in row: {!r}".format(row))
    if position < 1:
        raise ValueError("Position must be 1-based and positive: {!r}".format(row))
    return (str(row["chrom"]), position, str(row["ref"]).upper(), str(row["alt"]).upper())


def allele_key_string(row):
    chrom, position, ref, alt = allele_key(row)
    return "{}:{}:{}:{}".format(chrom, position, ref, alt)


def parse_bool(value, field="boolean"):
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    raise ValueError("Invalid {} value: {!r}".format(field, value))


def parse_optional_float(value, field="value"):
    if value is None or str(value).strip() in MISSING:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError("Invalid numeric {}: {!r}".format(field, value))
    if not math.isfinite(number):
        raise ValueError("Non-finite numeric {}: {!r}".format(field, value))
    return number


def corrected_from_medians(treated_median, control_median):
    """Return the frozen scientific target; missing group data stays missing."""
    if treated_median is None or control_median is None:
        return None
    return max(treated_median - control_median, 0.0)


def require_columns(fields, required, label):
    missing = [field for field in required if field not in fields]
    if missing:
        raise ValueError("{} is missing required columns: {}".format(label, ", ".join(missing)))


def index_unique(rows, label):
    result = {}
    for row in rows:
        key = allele_key(row)
        if key in result:
            raise ValueError("Duplicate allele key in {}: {}".format(label, key))
        result[key] = row
    return result


def close_enough(left, right, tolerance=1e-8):
    return math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance)


def validate_label_row(row):
    key_text = allele_key_string(row)
    treated = parse_optional_float(row.get("treated_median"), "treated_median")
    control = parse_optional_float(row.get("control_median"), "control_median")
    raw = parse_optional_float(row.get("raw_edit_rate_difference"), "raw_edit_rate_difference")
    corrected = parse_optional_float(
        row.get("corrected_editing_efficiency"), "corrected_editing_efficiency"
    )
    eligible = parse_bool(row.get("training_eligible"), "training_eligible")

    for name, value in (("treated_median", treated), ("control_median", control), ("corrected", corrected)):
        if value is not None and not 0.0 <= value <= 1.0:
            raise ValueError("{} outside [0,1] for {}: {}".format(name, key_text, value))

    if raw is None:
        if corrected is not None:
            raise ValueError("Missing raw difference has a non-missing corrected label for {}".format(key_text))
    else:
        if treated is None or control is None:
            raise ValueError("Raw difference requires both replicate medians for {}".format(key_text))
        expected_raw = treated - control
        if not close_enough(raw, expected_raw):
            raise ValueError(
                "Raw difference does not equal treated_median-control_median for {}: {} vs {}".format(
                    key_text, raw, expected_raw
                )
            )
        expected_corrected = max(raw, 0.0)
        if corrected is None or not close_enough(corrected, expected_corrected):
            raise ValueError(
                "Corrected label does not equal max(raw,0) for {}: {} vs {}".format(
                    key_text, corrected, expected_corrected
                )
            )

    if eligible and corrected is None:
        raise ValueError("Training-eligible row has a missing corrected label: {}".format(key_text))
    if not eligible:
        reason = str(row.get("exclusion_reason", "")).strip()
        if reason in MISSING or reason.lower() == "none":
            raise ValueError("Ineligible row lacks an explicit exclusion reason: {}".format(key_text))


def validate_sequence_row(row):
    key_text = allele_key_string(row)
    sequence = str(row.get("sequence_context", "")).upper()
    try:
        sequence_length = int(row.get("sequence_length", ""))
        center_index = int(row.get("center_index", ""))
    except (TypeError, ValueError):
        raise ValueError("Invalid sequence length/center index for {}".format(key_text))
    if sequence_length != SEQUENCE_LENGTH or len(sequence) != SEQUENCE_LENGTH:
        raise ValueError("Sequence must be exactly 101 nt for {}".format(key_text))
    if center_index != CENTER_INDEX:
        raise ValueError("Center index must be zero-based index 50 for {}".format(key_text))
    if sequence[CENTER_INDEX] != "C" or str(row.get("center_base", "")).upper() != "C":
        raise ValueError("Transcript-oriented center base is not C for {}".format(key_text))
    if str(row.get("transcript_oriented_ref", "")).upper() != "C":
        raise ValueError("Transcript-oriented REF is not C for {}".format(key_text))
    if str(row.get("transcript_oriented_alt", "")).upper() != "T":
        raise ValueError("Transcript-oriented ALT is not T for {}".format(key_text))
    if str(row.get("orientation_qc", "")).lower() != "pass":
        raise ValueError("Sequence orientation QC did not pass for {}".format(key_text))


def validate_shared_fields(label_row, metadata_row):
    key_text = allele_key_string(label_row)
    numeric = ("corrected_editing_efficiency", "raw_edit_rate_difference")
    for field in numeric:
        left = parse_optional_float(label_row.get(field), field)
        right = parse_optional_float(metadata_row.get(field), field)
        if left is None or right is None:
            if left is not None or right is not None:
                raise ValueError("Label/metadata missingness differs for {} at {}".format(field, key_text))
        elif not close_enough(left, right):
            raise ValueError("Label/metadata {} differs for {}".format(field, key_text))
    if parse_bool(label_row.get("training_eligible"), "training_eligible") != parse_bool(
        metadata_row.get("training_eligible"), "training_eligible"
    ):
        raise ValueError("Label/metadata eligibility differs for {}".format(key_text))
    for field in ("label_confidence", "exclusion_reason"):
        if str(label_row.get(field, "")) != str(metadata_row.get(field, "")):
            raise ValueError("Label/metadata {} differs for {}".format(field, key_text))


def merge_rows(label_fields, metadata_fields, label_by_key, metadata_by_key, construct=None):
    if set(label_by_key) != set(metadata_by_key):
        labels_only = sorted(set(label_by_key) - set(metadata_by_key))[:5]
        metadata_only = sorted(set(metadata_by_key) - set(label_by_key))[:5]
        raise ValueError(
            "Label and metadata keys are not an exact one-to-one match; labels-only={!r}, metadata-only={!r}".format(
                labels_only, metadata_only
            )
        )
    if construct is not None and set(construct) != set(label_by_key):
        raise ValueError("Construct metadata must contain exactly one row for every audited allele key")

    joined = []
    for key in sorted(label_by_key):
        label_row = label_by_key[key]
        metadata_row = metadata_by_key[key]
        validate_label_row(label_row)
        validate_sequence_row(metadata_row)
        validate_shared_fields(label_row, metadata_row)
        row = {
            "chrom": key[0],
            "position": key[1],
            "ref": key[2],
            "alt": key[3],
            "sequence": str(metadata_row["sequence_context"]).upper(),
        }
        for field in metadata_fields:
            if field not in row and field not in KEY_FIELDS:
                row[field] = metadata_row.get(field)
        for field in label_fields:
            if field not in row and field not in KEY_FIELDS:
                row[field] = label_row.get(field)
        if construct is not None:
            for field, value in construct[key].items():
                if field in KEY_FIELDS:
                    continue
                output_field = field if field not in row else "construct_{}".format(field)
                row[output_field] = value
        joined.append(row)
    return joined


def construct_subsets(rows):
    all_eligible = [row for row in rows if parse_bool(row["training_eligible"], "training_eligible")]
    excluded = [row for row in rows if not parse_bool(row["training_eligible"], "training_eligible")]
    high = [row for row in all_eligible if str(row.get("label_confidence", "")).lower() == "high"]
    high_low_control = [
        row
        for row in high
        if not parse_bool(row.get("elevated_control_background"), "elevated_control_background")
    ]
    if any(parse_optional_float(row.get("corrected_editing_efficiency")) is None for row in all_eligible):
        raise ValueError("Eligible subset contains a missing corrected label")
    if any(parse_bool(row["training_eligible"]) for row in excluded):
        raise AssertionError("Eligible row leaked into excluded output")
    if any(not parse_bool(row["training_eligible"]) for row in all_eligible + high + high_low_control):
        raise AssertionError("Ineligible row leaked into an eligible output")
    return {
        "all_eligible": all_eligible,
        "high_confidence": high,
        "high_confidence_low_control": high_low_control,
        "excluded": excluded,
    }


def build_overlap_clusters(rows):
    """Return allele-key -> merged 101-nt genomic-overlap cluster ID."""
    ordered = sorted(rows, key=lambda row: (str(row["chrom"]), int(row["position"]), str(row["ref"]), str(row["alt"])))
    clusters = {}
    cluster_number = 0
    current_chrom = None
    current_end = None
    current_id = None
    for row in ordered:
        chrom = str(row["chrom"])
        position = int(row["position"])
        start = position - CENTER_INDEX
        end = position + CENTER_INDEX
        if start < 1:
            raise ValueError("A valid 101-nt genomic window cannot start before position 1")
        if chrom != current_chrom or current_end is None or start > current_end:
            cluster_number += 1
            current_id = "overlap_cluster_{:06d}".format(cluster_number)
            current_chrom = chrom
            current_end = end
        else:
            current_end = max(current_end, end)
        clusters[allele_key(row)] = current_id
    return clusters


def build_duplicate_sequence_groups(rows):
    sequences = sorted({str(row["sequence"]).upper() for row in rows})
    id_by_sequence = {
        sequence: "duplicate_sequence_{:06d}".format(index + 1)
        for index, sequence in enumerate(sequences)
    }
    return {allele_key(row): id_by_sequence[str(row["sequence"]).upper()] for row in rows}


def build_gene_groups(rows):
    """Return allele-key -> unambiguous gene ID for gene-disjoint splitting."""
    groups = {}
    for row in rows:
        gene_id = str(row.get("gene_id", "")).strip()
        if gene_id in MISSING or "|" in gene_id:
            raise ValueError(
                "gene_disjoint splitting requires one unambiguous gene_id per row: {}".format(
                    allele_key_string(row)
                )
            )
        groups[allele_key(row)] = gene_id
    return groups


class UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, value):
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def union_groups(union_find, groups):
    for indices in groups.values():
        first = indices[0]
        for index in indices[1:]:
            union_find.union(first, index)


def label_class(row):
    target = parse_optional_float(row.get("corrected_editing_efficiency"), "corrected label")
    if target is None:
        raise ValueError("Split input contains a missing corrected label")
    return "zero" if close_enough(target, 0.0) else "positive"


def component_signature(component, rows):
    return "|".join(sorted(allele_key_string(rows[index]) for index in component))


def balanced_component_assignment(components, rows, seed):
    if len(components) < 3:
        raise ValueError("At least three independent split components are required")
    totals = {
        "rows": len(rows),
        "zero": sum(label_class(row) == "zero" for row in rows),
        "positive": sum(label_class(row) == "positive" for row in rows),
    }
    if totals["zero"] == 0 or totals["positive"] == 0:
        raise ValueError("Eligible data must contain both zero and positive corrected labels")
    targets = {
        split: {metric: SPLIT_FRACTIONS[split] * total for metric, total in totals.items()}
        for split in SPLIT_FRACTIONS
    }
    counts = {split: {"rows": 0, "zero": 0, "positive": 0} for split in SPLIT_FRACTIONS}
    assignments = {}

    component_stats = []
    for component in components:
        zero = sum(label_class(rows[index]) == "zero" for index in component)
        positive = len(component) - zero
        signature = component_signature(component, rows)
        tie = int(hashlib.sha256((str(seed) + "|" + signature).encode("utf-8")).hexdigest(), 16)
        component_stats.append((component, {"rows": len(component), "zero": zero, "positive": positive}, tie))
    component_stats.sort(key=lambda item: (-item[1]["rows"], -max(item[1]["zero"], item[1]["positive"]), item[2]))

    split_order = list(SPLIT_FRACTIONS)
    random.Random(seed).shuffle(split_order)

    for component, stats, tie in component_stats:
        candidates = []
        for split_index, split in enumerate(split_order):
            score = 0.0
            for candidate_split in SPLIT_FRACTIONS:
                for metric in ("rows", "zero", "positive"):
                    observed = counts[candidate_split][metric]
                    if candidate_split == split:
                        observed += stats[metric]
                    target = targets[candidate_split][metric]
                    score += ((observed - target) / max(target, 1.0)) ** 2
                    if observed > target * 1.20 and target >= 1:
                        score += 0.25 * ((observed - target) / target) ** 2
            deterministic_jitter = ((tie >> (split_index * 8)) & 0xFF) * 1e-12
            candidates.append((score + deterministic_jitter, split))
        chosen = min(candidates)[1]
        assignments[id(component)] = chosen
        for metric in counts[chosen]:
            counts[chosen][metric] += stats[metric]

    # With small synthetic inputs, protect the required train class coverage by
    # moving the smallest safe component. Production inputs normally need no move.
    for required_class in ("zero", "positive"):
        if counts["train"][required_class] > 0:
            continue
        movable = []
        for component, stats, tie in component_stats:
            source = assignments[id(component)]
            if source == "train" or stats[required_class] == 0:
                continue
            if counts[source]["rows"] - stats["rows"] <= 0:
                continue
            movable.append((stats["rows"], tie, component, stats, source))
        if not movable:
            raise ValueError("Training split cannot be given a {} label without emptying another split".format(required_class))
        _, _, component, stats, source = min(movable)
        assignments[id(component)] = "train"
        for metric in counts[source]:
            counts[source][metric] -= stats[metric]
            counts["train"][metric] += stats[metric]

    if any(counts[split]["rows"] == 0 for split in SPLIT_FRACTIONS):
        raise ValueError("Deterministic component allocation produced an empty split")
    return assignments, counts, component_stats


def assign_splits(rows, strategy="overlap_cluster", seed=DEFAULT_SEED):
    if strategy not in {"overlap_cluster", "chromosome", "gene_disjoint"}:
        raise ValueError("Unknown split strategy: {}".format(strategy))
    rows = sorted(rows, key=lambda row: allele_key(row))
    if len(index_unique(rows, "eligible split rows")) != len(rows):
        raise AssertionError("Duplicate split keys")
    overlap_by_key = build_overlap_clusters(rows)
    sequence_group_by_key = build_duplicate_sequence_groups(rows)
    gene_group_by_key = build_gene_groups(rows) if strategy == "gene_disjoint" else {}
    union_find = UnionFind(len(rows))

    primary_groups = defaultdict(list)
    sequence_groups = defaultdict(list)
    for index, row in enumerate(rows):
        key = allele_key(row)
        if strategy == "overlap_cluster":
            primary_id = overlap_by_key[key]
        elif strategy == "chromosome":
            primary_id = str(row["chrom"])
        else:
            primary_id = gene_group_by_key[key]
        primary_groups[primary_id].append(index)
        sequence_groups[sequence_group_by_key[key]].append(index)
    union_groups(union_find, primary_groups)
    union_groups(union_find, sequence_groups)
    if strategy == "gene_disjoint":
        overlap_groups = defaultdict(list)
        for index, row in enumerate(rows):
            overlap_groups[overlap_by_key[allele_key(row)]].append(index)
        union_groups(union_find, overlap_groups)

    components_by_root = defaultdict(list)
    for index in range(len(rows)):
        components_by_root[union_find.find(index)].append(index)
    components = list(components_by_root.values())
    assignments, _, component_stats = balanced_component_assignment(components, rows, seed)

    split_rows = []
    for component, _, _ in component_stats:
        split = assignments[id(component)]
        for index in component:
            row = rows[index]
            key = allele_key(row)
            split_rows.append(
                {
                    "allele_key": allele_key_string(row),
                    "chrom": key[0],
                    "position": key[1],
                    "ref": key[2],
                    "alt": key[3],
                    "sequence": row["sequence"],
                    "center_index": CENTER_INDEX,
                    "corrected_editing_efficiency": row["corrected_editing_efficiency"],
                    "split": split,
                    "overlap_cluster_id": overlap_by_key[key],
                    "duplicate_sequence_group_id": sequence_group_by_key[key],
                    "gene_group_id": gene_group_by_key.get(key, str(row.get("gene_id", "NA"))),
                    "label_confidence": row["label_confidence"],
                    "label_class": label_class(row),
                    "training_eligible": row["training_eligible"],
                    "elevated_control_background": row["elevated_control_background"],
                    "raw_edit_rate_difference": row["raw_edit_rate_difference"],
                    "control_median": row["control_median"],
                    "treated_median": row["treated_median"],
                }
            )
    split_rows.sort(key=lambda row: (row["split"], row["chrom"], int(row["position"]), row["ref"], row["alt"]))
    return split_rows


def split_leakage_count(rows, group_field):
    splits_by_group = defaultdict(set)
    for row in rows:
        splits_by_group[str(row[group_field])].add(str(row["split"]))
    return sum(len(splits) > 1 for splits in splits_by_group.values())


def validate_split_assignments(rows, strategy="overlap_cluster"):
    keys = [allele_key(row) for row in rows]
    duplicate_key_count = len(keys) - len(set(keys))
    overlap_leakage = split_leakage_count(rows, "overlap_cluster_id")
    sequence_leakage = split_leakage_count(rows, "duplicate_sequence_group_id")
    split_names = {str(row["split"]) for row in rows}
    missing_splits = set(SPLIT_FRACTIONS) - split_names
    train = [row for row in rows if row["split"] == "train"]
    train_classes = {row["label_class"] for row in train}
    chromosome_leakage = 0
    gene_leakage = 0
    if strategy == "chromosome":
        splits_by_chrom = defaultdict(set)
        for row in rows:
            splits_by_chrom[str(row["chrom"])].add(str(row["split"]))
        chromosome_leakage = sum(len(splits) > 1 for splits in splits_by_chrom.values())
    if strategy == "gene_disjoint":
        gene_leakage = split_leakage_count(rows, "gene_group_id")
    errors = []
    if duplicate_key_count:
        errors.append("allele key appears more than once")
    if overlap_leakage:
        errors.append("overlap cluster crosses splits")
    if sequence_leakage:
        errors.append("identical sequence crosses splits")
    if chromosome_leakage:
        errors.append("chromosome crosses splits")
    if gene_leakage:
        errors.append("gene crosses splits")
    if missing_splits:
        errors.append("empty split(s): {}".format(",".join(sorted(missing_splits))))
    if "zero" not in train_classes or "positive" not in train_classes:
        errors.append("training split lacks zero or positive labels")
    if errors:
        raise ValueError("Invalid leakage-resistant split: " + "; ".join(errors))
    return {
        "allele_key_leakage_count": duplicate_key_count,
        "overlap_cluster_leakage_count": overlap_leakage,
        "duplicate_sequence_leakage_count": sequence_leakage,
        "chromosome_leakage_count": chromosome_leakage,
        "gene_leakage_count": gene_leakage,
    }


def subset_counts(rows, selector):
    selected = [row for row in rows if selector(row)]
    result = {}
    for split in SPLIT_FRACTIONS:
        group = [row for row in selected if row["split"] == split]
        result[split] = {
            "rows": len(group),
            "zero": sum(row["label_class"] == "zero" for row in group),
            "positive": sum(row["label_class"] == "positive" for row in group),
        }
    return result


def make_split_qc(split_rows, strategy, seed):
    leakage = validate_split_assignments(split_rows, strategy)
    counts = subset_counts(split_rows, lambda row: True)
    high_counts = subset_counts(
        split_rows, lambda row: str(row["label_confidence"]).lower() == "high"
    )
    if any(high_counts[split]["rows"] == 0 for split in SPLIT_FRACTIONS):
        raise ValueError("Recommended high-confidence dataset has an empty split")
    if high_counts["train"]["zero"] == 0 or high_counts["train"]["positive"] == 0:
        raise ValueError("Recommended high-confidence training split lacks zero or positive labels")
    return {
        "schema_version": 1,
        "seed": seed,
        "strategy": strategy,
        "target_fractions": SPLIT_FRACTIONS,
        "algorithm": (
            "Union overlapping 101-nt genomic windows with exact-sequence groups, then assign connected components "
            "by deterministic multi-objective greedy balancing"
            if strategy == "overlap_cluster"
            else (
                "Union genes, overlapping 101-nt windows and exact-sequence groups, then assign connected components by deterministic multi-objective greedy balancing"
                if strategy == "gene_disjoint"
                else "Union chromosomes with exact-sequence groups, then assign connected components by deterministic multi-objective greedy balancing"
            )
        ),
        "all_eligible_counts": counts,
        "high_confidence_counts": high_counts,
        "overlap_cluster_count": len({row["overlap_cluster_id"] for row in split_rows}),
        "duplicate_sequence_group_count": len({row["duplicate_sequence_group_id"] for row in split_rows}),
        "gene_group_count": len({row["gene_group_id"] for row in split_rows}),
        "unique_sequence_count": len({row["sequence"] for row in split_rows}),
        "leakage_checks": leakage,
        "validation_status": "pass",
    }


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_checksums(directory, filenames=None):
    directory = Path(directory)
    if filenames is None:
        paths = sorted(
            path for path in directory.iterdir() if path.is_file() and path.name != "checksums.sha256"
        )
    else:
        paths = [directory / name for name in filenames]
    checksum_path = directory / "checksums.sha256"
    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in paths:
            if not path.is_file():
                raise FileNotFoundError("Cannot checksum missing file: {}".format(path))
            handle.write("{}  {}\n".format(sha256_file(path), path.name))
    return checksum_path


def verify_checksums(directory):
    directory = Path(directory)
    checksum_path = directory / "checksums.sha256"
    with checksum_path.open(encoding="utf-8") as handle:
        for line in handle:
            expected, filename = line.rstrip("\n").split("  ", 1)
            observed = sha256_file(directory / filename)
            if observed != expected:
                raise ValueError("Checksum mismatch for {}".format(filename))
    return True


def field_type(field):
    if (
        field
        in {
            "position",
            "sequence_length",
            "center_index",
            "minimum_usable_depth",
            "minimum_covered_replicates_per_group",
        }
        or field.endswith("_count")
        or field.endswith("_replicates")
        or field.endswith("_usable_depth")
    ):
        return "integer"
    if field in {
        "corrected_editing_efficiency",
        "raw_edit_rate_difference",
        "control_median",
        "treated_median",
        "treated_mean",
        "control_mean",
        "treated_MAD",
        "control_MAD",
        "treated_range",
        "control_range",
        "BH_FDR",
        "fisher_exact_screening_pvalue",
    } or field.endswith("_edit_rate"):
        return "number_or_NA"
    if field in {
        "training_eligible",
        "elevated_control_background",
        "passes_FDR_0.05",
        "treated_replicate_consistent",
        "control_replicate_consistent",
        "sufficient_coverage_all_six",
    }:
        return "boolean_0_or_1"
    return "string"


FIELD_DESCRIPTIONS = {
    "allele_key": "Exact GRCh38 genomic allele key chrom:position:ref:alt.",
    "chrom": "GRCh38 chromosome or contig name.",
    "position": "GRCh38 1-based edited genomic position.",
    "ref": "Genomic reference allele.",
    "alt": "Genomic alternate allele.",
    "sequence": "101-nt transcript-oriented sequence centered on the edited cytidine.",
    "sequence_context": "Source 101-nt transcript-oriented sequence; identical to sequence.",
    "corrected_editing_efficiency": "max(treated replicate-rate median - control replicate-rate median, 0).",
    "raw_edit_rate_difference": "Treated replicate-rate median minus control replicate-rate median before clipping.",
    "split": "Leakage-resistant train, validation, or test assignment.",
    "overlap_cluster_id": "Merged same-chromosome cluster of overlapping position±50 intervals.",
    "duplicate_sequence_group_id": "Group identifier shared by exactly identical 101-nt sequences.",
    "gene_group_id": "Unambiguous gene identifier used for gene-disjoint assignment when requested.",
    "label_class": "zero if corrected target is 0; positive otherwise.",
    "training_eligible": "Frozen audit flag; requires at least 2/3 sufficiently covered replicates in each group plus audit QC.",
    "label_confidence": "Frozen audit confidence; high requires 3+3 covered replicates and replicate-consistency thresholds.",
    "elevated_control_background": "Frozen audit flag for control median at or above the documented threshold.",
    "exclusion_reason": "Explicit semicolon-delimited reason an audited row is not training eligible.",
    "center_index": "Zero-based nucleotide index 50 in the un-tokenized 101-nt sequence.",
}


def write_data_dictionary(path, joined_fields, split_fields):
    rows = []
    for table, fields in (
        ("eligible_and_excluded_tables", joined_fields),
        ("CU5.17_lamar_splits.tsv.gz", split_fields),
    ):
        for field in fields:
            rows.append(
                {
                    "table": table,
                    "column": field,
                    "type": field_type(field),
                    "description": FIELD_DESCRIPTIONS.get(
                        field, "Audited source field retained for reproducibility and label traceability."
                    ),
                }
            )
    write_tsv(path, ("table", "column", "type", "description"), rows)


def write_readme(path, counts, strategy, seed):
    text = """# CU5.17 LAMAR fine-tuning handoff

## Recommended primary analysis

Use `CU5.17_lamar_high_confidence.tsv.gz` as the primary biological dataset and
use the high-confidence rows in `CU5.17_lamar_splits.tsv.gz` for model fitting.
High confidence retains the frozen requirement for all three treated and all
three control replicates plus the documented replicate-consistency thresholds.

## Sensitivity analyses

- `CU5.17_lamar_all_eligible.tsv.gz` retains every row passing the frozen audit
  (at least 2/3 sufficiently covered replicates in each group).
- `CU5.17_lamar_high_confidence_low_control.tsv.gz` removes high-confidence rows
  flagged for elevated control background and is a stricter sensitivity set.
- `CU5.17_lamar_excluded.tsv.gz` is audit-only and preserves explicit reasons;
  it must not be silently relabeled or included as eligible data.

The split table covers all eligible rows so primary and sensitivity analyses use
the same leakage-resistant partition. Filter `label_confidence == high` for the
primary analysis. Zero-valued corrected targets are retained as valid examples.

## Frozen target and data boundary

`corrected_editing_efficiency = max(treated_median - control_median, 0)`.
Missing control evidence is never converted to zero. The 3,333 final candidates
are a subset of the broad 9,930-site universe and are not an independent test
set. The broad universe is itself derived from called candidate sites; it is not
a complete transcriptome-wide negative universe.

This handoff does not contain or invent `puf_target_seq`, `label_total_count`, or
non-center token labels. Computational QC is not experimental validation.

## Split

- strategy: `{strategy}`
- deterministic seed: `{seed}`
- target allocation: 80/10/10 train/validation/test
- all eligible: `{all_eligible}` rows
- high confidence: `{high_confidence}` rows
- high confidence with low control background: `{high_confidence_low_control}` rows
- excluded: `{excluded}` rows

`overlap_cluster` merges genomic position±50 intervals on each chromosome and
also links identical sequences across coordinates before component-level split
assignment. `chromosome` instead holds out whole chromosomes (while still
linking identical sequences) and therefore measures a stronger genomic
distribution shift.

`gene_disjoint` keeps every gene in exactly one split and also links overlapping
windows and identical sequences. Rows with missing or ambiguous `gene_id` stop
the build instead of being assigned by coordinate.

## Scalar export

```bash
python pipeline/scripts/rna/export_lamar_scalar_regression.py \\
  --input CU5.17_lamar_splits.tsv.gz \\
  --output CU5.17_lamar_scalar_high_confidence.tsv.gz
```

The exported center index is a nucleotide index. Confirm how the selected LAMAR
tokenizer adds special tokens before mapping it to a model-token index.
""".format(strategy=strategy, seed=seed, **counts)
    Path(path).write_text(text, encoding="utf-8")


def safe_output_directory(output_dir, overwrite=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    existing = [output_dir / name for name in EXPECTED_OUTPUTS if (output_dir / name).exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Output files already exist; use --overwrite to replace only known handoff files: {}".format(
                ", ".join(path.name for path in existing)
            )
        )
    if overwrite:
        for path in existing:
            if path.is_file() or path.is_symlink():
                path.unlink()
            else:
                raise ValueError("Refusing to replace non-file output: {}".format(path))
    return output_dir


def copy_public_outputs(source, destination, overwrite=False):
    source = Path(source)
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    for name in PUBLIC_OUTPUTS:
        target = destination / name
        if target.exists() and not overwrite:
            raise FileExistsError("Public output already exists: {}".format(target))
        shutil.copyfile(source / name, target)
    write_checksums(destination, PUBLIC_OUTPUTS)
    verify_checksums(destination)


def build_handoff(
    labels_path,
    metadata_path,
    output_dir,
    construct_metadata_path=None,
    seed=DEFAULT_SEED,
    split_strategy="overlap_cluster",
    overwrite=False,
    public_copy_dir=None,
):
    output_dir = safe_output_directory(output_dir, overwrite=overwrite)
    label_fields, label_rows = read_tsv(labels_path)
    metadata_fields, metadata_rows = read_tsv(metadata_path)
    label_required = KEY_FIELDS + (
        "treated_median",
        "control_median",
        "raw_edit_rate_difference",
        "corrected_editing_efficiency",
        "training_eligible",
        "label_confidence",
        "elevated_control_background",
        "exclusion_reason",
    )
    metadata_required = KEY_FIELDS + (
        "sequence_context",
        "sequence_length",
        "center_index",
        "center_base",
        "orientation_qc",
        "transcript_oriented_ref",
        "transcript_oriented_alt",
        "corrected_editing_efficiency",
        "raw_edit_rate_difference",
        "training_eligible",
        "label_confidence",
        "exclusion_reason",
    )
    require_columns(label_fields, label_required, "background-corrected labels")
    require_columns(metadata_fields, metadata_required, "LAMAR-ready metadata")
    label_by_key = index_unique(label_rows, "background-corrected labels")
    metadata_by_key = index_unique(metadata_rows, "LAMAR-ready metadata")

    construct_by_key = None
    construct_name = None
    construct_sha = None
    if construct_metadata_path is not None:
        construct_fields, construct_rows = read_tsv(construct_metadata_path)
        require_columns(construct_fields, KEY_FIELDS, "construct metadata")
        construct_by_key = index_unique(construct_rows, "construct metadata")
        construct_name = Path(construct_metadata_path).name
        construct_sha = sha256_file(construct_metadata_path)

    joined = merge_rows(
        label_fields,
        metadata_fields,
        label_by_key,
        metadata_by_key,
        construct=construct_by_key,
    )
    subsets = construct_subsets(joined)
    if not subsets["high_confidence"]:
        raise ValueError("No high-confidence rows are available for the recommended primary dataset")

    joined_fields = list(joined[0])
    write_tsv(output_dir / "CU5.17_lamar_all_eligible.tsv.gz", joined_fields, subsets["all_eligible"])
    write_tsv(output_dir / "CU5.17_lamar_high_confidence.tsv.gz", joined_fields, subsets["high_confidence"])
    write_tsv(
        output_dir / "CU5.17_lamar_high_confidence_low_control.tsv.gz",
        joined_fields,
        subsets["high_confidence_low_control"],
    )
    write_tsv(output_dir / "CU5.17_lamar_excluded.tsv.gz", joined_fields, subsets["excluded"])

    split_rows = assign_splits(subsets["all_eligible"], split_strategy, seed)
    split_fields = list(split_rows[0])
    write_tsv(output_dir / "CU5.17_lamar_splits.tsv.gz", split_fields, split_rows)
    split_qc = make_split_qc(split_rows, split_strategy, seed)
    (output_dir / "split_qc.json").write_text(
        json.dumps(split_qc, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_data_dictionary(output_dir / "data_dictionary.tsv", joined_fields, split_fields)

    counts = {name: len(rows) for name, rows in subsets.items()}
    manifest = {
        "schema_version": 1,
        "builder": "pipeline/scripts/rna/prepare_lamar_finetuning_handoff.py",
        "packaging_runtime": {
            "python_version": sys.version.split()[0],
            "scientific_counts_recomputed": False,
            "bam_or_pileup_accessed": False,
        },
        "scientific_target": "max(median(treated editing rates) - median(control editing rates), 0)",
        "seed": seed,
        "split_strategy": split_strategy,
        "split_algorithm": split_qc["algorithm"],
        "split_fractions": SPLIT_FRACTIONS,
        "input_files": {
            "background_corrected_labels": {
                "filename": Path(labels_path).name,
                "sha256": sha256_file(labels_path),
                "rows": len(label_rows),
            },
            "lamar_ready_metadata": {
                "filename": Path(metadata_path).name,
                "sha256": sha256_file(metadata_path),
                "rows": len(metadata_rows),
            },
            "construct_metadata": (
                {"filename": construct_name, "sha256": construct_sha, "rows": len(construct_by_key)}
                if construct_by_key is not None
                else None
            ),
        },
        "row_counts": counts,
        "split_qc_file": "split_qc.json",
        "split_qc": split_qc,
        "recommended_primary": "CU5.17_lamar_high_confidence.tsv.gz",
        "sensitivity_datasets": [
            "CU5.17_lamar_all_eligible.tsv.gz",
            "CU5.17_lamar_high_confidence_low_control.tsv.gz",
        ],
        "split_population": "all training-eligible rows; primary analysis filters label_confidence == high",
        "guardrails": {
            "missing_labels_are_not_zero": True,
            "zero_corrected_labels_retained": True,
            "puf_target_seq_fabricated": False,
            "label_total_count_fabricated": False,
            "non_center_token_labels_created": False,
            "final_3333_used_as_independent_test": False,
        },
        "validation": {
            "exact_one_to_one_label_metadata_join": True,
            "sequence_length_101_center_index_50_center_C": True,
            "clipping_formula_checked": True,
            "leakage_checks": split_qc["leakage_checks"],
        },
    }
    (output_dir / "handoff_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_readme(output_dir / "README.md", counts, split_strategy, seed)
    write_checksums(output_dir)
    verify_checksums(output_dir)

    if public_copy_dir is not None:
        copy_public_outputs(output_dir, public_copy_dir, overwrite=overwrite)
    return {"counts": counts, "split_qc": split_qc, "manifest": manifest}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, required=True, help="Frozen background_corrected_labels.tsv.gz")
    parser.add_argument("--metadata", type=Path, required=True, help="Frozen lamar_ready_metadata.tsv.gz")
    parser.add_argument("--construct-metadata", type=Path, help="Optional exact one-row-per-allele construct metadata")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--split-strategy", choices=("overlap_cluster", "chromosome", "gene_disjoint"), default="overlap_cluster"
    )
    parser.add_argument(
        "--public-copy-dir",
        type=Path,
        help="Optionally copy only compact public model-facing outputs and their checksums",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace only the known handoff output files")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        result = build_handoff(
            labels_path=args.labels,
            metadata_path=args.metadata,
            construct_metadata_path=args.construct_metadata,
            output_dir=args.output_dir,
            seed=args.seed,
            split_strategy=args.split_strategy,
            overwrite=args.overwrite,
            public_copy_dir=args.public_copy_dir,
        )
    except Exception as error:
        raise SystemExit("ERROR: {}".format(error))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
