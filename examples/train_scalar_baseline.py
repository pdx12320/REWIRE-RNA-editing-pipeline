#!/usr/bin/env python3
"""Run leakage-safe scalar baselines on the precomputed CU5.17 split table.

This is a data-integrity check, not a final biological model.  The script never
creates a new random split.  By default it evaluates the recommended
high-confidence subset using the assignments already present in the handoff.
"""

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path


SPLITS = ("train", "validation", "test")


def open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def parse_bool(value):
    text = str(value).lower()
    if text in {"1", "true", "yes"}:
        return True
    if text in {"0", "false", "no"}:
        return False
    raise ValueError("Invalid boolean value: {!r}".format(value))


def read_split_rows(path, subset):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "chrom",
            "position",
            "ref",
            "alt",
            "sequence",
            "corrected_editing_efficiency",
            "split",
            "overlap_cluster_id",
            "duplicate_sequence_group_id",
            "label_confidence",
            "training_eligible",
            "elevated_control_background",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Split file is missing: {}".format(", ".join(sorted(missing))))
        source_rows = [dict(row) for row in reader]

    seen = set()
    cluster_splits = defaultdict(set)
    sequence_splits = defaultdict(set)
    rows = []
    for row in source_rows:
        key = (row["chrom"], int(row["position"]), row["ref"], row["alt"])
        if key in seen:
            raise ValueError("Duplicate allele key: {}".format(key))
        seen.add(key)
        cluster_splits[row["overlap_cluster_id"]].add(row["split"])
        sequence_splits[row["duplicate_sequence_group_id"]].add(row["split"])
        if not parse_bool(row["training_eligible"]):
            raise ValueError("Ineligible row found in split file: {}".format(key))
        include = subset == "all_eligible"
        if subset in {"high_confidence", "high_confidence_low_control"}:
            include = row["label_confidence"].lower() == "high"
        if subset == "high_confidence_low_control" and include:
            include = not parse_bool(row["elevated_control_background"])
        if include:
            row["target"] = float(row["corrected_editing_efficiency"])
            rows.append(row)
    if any(len(splits) > 1 for splits in cluster_splits.values()):
        raise ValueError("Overlap-cluster leakage detected")
    if any(len(splits) > 1 for splits in sequence_splits.values()):
        raise ValueError("Duplicate-sequence leakage detected")
    if {row["split"] for row in rows} != set(SPLITS):
        raise ValueError("Selected subset must have non-empty train/validation/test splits")
    train_targets = [row["target"] for row in rows if row["split"] == "train"]
    if not any(value == 0 for value in train_targets) or not any(value > 0 for value in train_targets):
        raise ValueError("Training split must include zero and positive labels")
    return rows


def average_ranks(values):
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][1] == ordered[index][1]:
            end += 1
        average = (index + 1 + end) / 2.0
        for cursor in range(index, end):
            ranks[ordered[cursor][0]] = average
        index = end
    return ranks


def pearson(left, right):
    if len(left) < 2:
        return None
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_norm = math.sqrt(sum((a - left_mean) ** 2 for a in left))
    right_norm = math.sqrt(sum((b - right_mean) ** 2 for b in right))
    if left_norm == 0 or right_norm == 0:
        return None
    return numerator / (left_norm * right_norm)


def metrics(observed, predicted):
    if not observed:
        return {"n": 0, "mae": None, "rmse": None, "pearson": None, "spearman": None}
    errors = [prediction - truth for truth, prediction in zip(observed, predicted)]
    return {
        "n": len(observed),
        "mae": sum(abs(error) for error in errors) / len(errors),
        "rmse": math.sqrt(sum(error * error for error in errors) / len(errors)),
        "pearson": pearson(observed, predicted),
        "spearman": pearson(average_ranks(observed), average_ranks(predicted)),
    }


def stratified_metrics(rows, predictions):
    result = {}
    scopes = {
        "all": list(range(len(rows))),
        "zero_label": [index for index, row in enumerate(rows) if row["target"] == 0],
        "positive_label": [index for index, row in enumerate(rows) if row["target"] > 0],
    }
    for name, indices in scopes.items():
        result[name] = metrics(
            [rows[index]["target"] for index in indices],
            [predictions[index] for index in indices],
        )
    return result


def run_baselines(rows, kmer_size=3, alpha=1.0):
    try:
        from sklearn.feature_extraction.text import CountVectorizer
        from sklearn.linear_model import Ridge
    except ImportError:
        raise RuntimeError(
            "scikit-learn is required for the k-mer ridge baseline; create "
            "pipeline/env/lamar_scalar_baseline.yml"
        )

    by_split = {split: [row for row in rows if row["split"] == split] for split in SPLITS}
    train_mean = sum(row["target"] for row in by_split["train"]) / len(by_split["train"])
    vectorizer = CountVectorizer(
        analyzer="char", ngram_range=(kmer_size, kmer_size), lowercase=False, dtype=float
    )
    train_matrix = vectorizer.fit_transform([row["sequence"] for row in by_split["train"]])
    model = Ridge(alpha=alpha)
    model.fit(train_matrix, [row["target"] for row in by_split["train"]])

    report = {
        "purpose": "data-integrity baseline; not a final biological model",
        "split_source": "precomputed input assignments; no random row split performed",
        "kmer_size": kmer_size,
        "ridge_alpha": alpha,
        "row_counts": {split: len(by_split[split]) for split in SPLITS},
        "metrics": {"constant_train_mean": {}, "kmer_ridge": {}},
    }
    for split in SPLITS:
        split_rows = by_split[split]
        constant_predictions = [train_mean] * len(split_rows)
        matrix = vectorizer.transform([row["sequence"] for row in split_rows])
        ridge_predictions = [min(max(float(value), 0.0), 1.0) for value in model.predict(matrix)]
        report["metrics"]["constant_train_mean"][split] = stratified_metrics(
            split_rows, constant_predictions
        )
        report["metrics"]["kmer_ridge"][split] = stratified_metrics(split_rows, ridge_predictions)
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="CU5.17_lamar_splits.tsv.gz")
    parser.add_argument(
        "--subset",
        choices=("high_confidence", "all_eligible", "high_confidence_low_control"),
        default="high_confidence",
    )
    parser.add_argument("--kmer-size", type=int, default=3)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--output-json", type=Path)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        rows = read_split_rows(args.input, args.subset)
        report = run_baselines(rows, args.kmer_size, args.alpha)
    except Exception as error:
        raise SystemExit("ERROR: {}".format(error))
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.output_json:
        args.output_json.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
