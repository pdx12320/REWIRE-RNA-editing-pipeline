#!/usr/bin/env python3
"""Filter mappability-audited C sites, enforce 1:N balance and create gene-disjoint splits."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import prepare_lamar_finetuning_handoff as handoff


def open_text(path: Path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_rows(path: Path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Mappability audit table has no header")
        rows = [dict(row) for row in reader]
    required = {
        "chrom", "position", "ref", "alt", "gene_id", "sequence_context",
        "label_class", "corrected_editing_efficiency", "training_eligible", "mappability_pass",
    }
    missing = sorted(required - set(reader.fieldnames))
    if missing:
        raise ValueError("Mappability table is missing columns: {}".format(", ".join(missing)))
    return list(reader.fieldnames), rows


def allele_key(row):
    return "{}:{}:{}:{}".format(row["chrom"], int(row["position"]), row["ref"], row["alt"])


def deterministic_rank(row, seed: int):
    return hashlib.sha256((str(seed) + "|" + allele_key(row)).encode("utf-8")).hexdigest()


def truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def select_rows(rows, negative_ratio: int, seed: int):
    passed = [row for row in rows if truthy(row["training_eligible"]) and truthy(row["mappability_pass"])]
    positives = [row for row in passed if row["label_class"] == "positive"]
    negatives = [row for row in passed if row["label_class"] == "strict_negative"]
    unexpected = sorted({row["label_class"] for row in passed} - {"positive", "strict_negative"})
    if unexpected:
        raise ValueError("Training-eligible rows contain unexpected classes: {}".format(", ".join(unexpected)))
    if not positives:
        raise ValueError("No mappability-passing positive sites remain")
    target_negatives = len(positives) * negative_ratio
    negatives.sort(key=lambda row: deterministic_rank(row, seed))
    selected_negatives = negatives[:target_negatives]
    selected = positives + selected_negatives
    for row in selected:
        row["training_eligible"] = "1"
        row["label_confidence"] = "high"
        row["exclusion_reason"] = "none"
    selected.sort(key=lambda row: (row["chrom"], int(row["position"]), row["ref"], row["alt"]))
    return selected, {
        "mappability_pass_positive": len(positives),
        "mappability_pass_strict_negative": len(negatives),
        "selected_positive": len(positives),
        "selected_strict_negative": len(selected_negatives),
        "requested_negative_ratio": negative_ratio,
        "realized_negative_ratio": len(selected_negatives) / len(positives),
        "negative_pool_exhausted": len(selected_negatives) < target_negatives,
    }


def finalize(input_path: Path, output_dir: Path, negative_ratio: int, seed: int, overwrite: bool):
    fields, rows = read_rows(input_path)
    selected, selection = select_rows(rows, negative_ratio, seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    labels = output_dir / "strict_selected_labels.tsv.gz"
    metadata = output_dir / "strict_selected_metadata.tsv.gz"
    for path in (labels, metadata):
        if path.exists() and not overwrite:
            raise FileExistsError("Output exists; use --overwrite: {}".format(path))
    handoff.write_tsv(labels, fields, selected)
    handoff.write_tsv(metadata, fields, selected)
    handoff_dir = output_dir / "gene_disjoint_handoff"
    result = handoff.build_handoff(
        labels,
        metadata,
        handoff_dir,
        seed=seed,
        split_strategy="gene_disjoint",
        overwrite=overwrite,
    )
    report = {
        "status": "pass",
        "input": str(input_path.resolve()),
        "selection": selection,
        "split_strategy": "gene_disjoint",
        "handoff": result,
    }
    report_path = output_dir / "strict_selection_audit.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mappability-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--negative-ratio", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260822)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    if args.negative_ratio < 1:
        parser.error("--negative-ratio must be positive")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        report = finalize(
            args.mappability_table,
            args.output_dir,
            args.negative_ratio,
            args.seed,
            args.overwrite,
        )
    except Exception as error:
        raise SystemExit("ERROR: {}".format(error))
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
