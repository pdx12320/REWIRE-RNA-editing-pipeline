#!/usr/bin/env python3
"""Export leakage-safe center-C scalar regression records for a LAMAR encoder.

The output is designed for a 101-nt sequence -> pretrained encoder -> center
hidden state -> scalar regression head.  ``center_index`` is a nucleotide index,
not a model-token index.  Tokenizers that prepend [CLS] or other special tokens
can shift the corresponding token index; callers must validate that mapping.

This is not a claim of native compatibility with the historical token-level
LAMAR trainer.
"""

import argparse
import csv
import gzip
import io
import json
import math
from pathlib import Path


REQUIRED_FIELDS = (
    "chrom",
    "position",
    "ref",
    "alt",
    "sequence",
    "center_index",
    "corrected_editing_efficiency",
    "split",
    "overlap_cluster_id",
    "duplicate_sequence_group_id",
    "label_confidence",
    "raw_edit_rate_difference",
    "control_median",
    "treated_median",
    "training_eligible",
    "elevated_control_background",
)

OUTPUT_FIELDS = (
    "sequence",
    "target",
    "center_index",
    "split",
    "overlap_cluster_id",
    "duplicate_sequence_group_id",
    "chrom",
    "position",
    "ref",
    "alt",
    "label_confidence",
    "raw_edit_rate_difference",
    "control_median",
    "treated_median",
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
            raise ValueError("Input has no TSV header: {}".format(path))
        missing = [field for field in REQUIRED_FIELDS if field not in reader.fieldnames]
        if missing:
            raise ValueError("Split table is missing: {}".format(", ".join(missing)))
        return [dict(row) for row in reader]


def write_tsv(path, fields, rows):
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
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        handle.close()
        if raw is not None:
            raw.close()


def parse_bool(value, field):
    text = str(value).strip().lower()
    if text in {"1", "true", "yes"}:
        return True
    if text in {"0", "false", "no"}:
        return False
    raise ValueError("Invalid {}: {!r}".format(field, value))


def select_rows(rows, subset):
    selected = []
    keys = set()
    for row in rows:
        key = (row["chrom"], int(row["position"]), row["ref"], row["alt"])
        if key in keys:
            raise ValueError("Duplicate allele key in split input: {}".format(key))
        keys.add(key)
        if not parse_bool(row["training_eligible"], "training_eligible"):
            raise ValueError("Split input contains an ineligible row: {}".format(key))
        include = subset == "all_eligible"
        if subset in {"high_confidence", "high_confidence_low_control"}:
            include = row["label_confidence"].lower() == "high"
        if subset == "high_confidence_low_control" and include:
            include = not parse_bool(
                row["elevated_control_background"], "elevated_control_background"
            )
        if include:
            sequence = row["sequence"].upper()
            center_index = int(row["center_index"])
            target = float(row["corrected_editing_efficiency"])
            if len(sequence) != 101 or center_index != 50 or sequence[center_index] != "C":
                raise ValueError("Invalid 101-nt center-C sequence for {}".format(key))
            if not math.isfinite(target) or not 0 <= target <= 1:
                raise ValueError("Invalid corrected target for {}".format(key))
            selected.append(row)
    if not selected:
        raise ValueError("Requested subset is empty: {}".format(subset))
    if {row["split"] for row in selected} != {"train", "validation", "test"}:
        raise ValueError("Requested subset must retain non-empty train/validation/test splits")
    return selected


def scalar_rows(rows):
    return [
        {
            "sequence": row["sequence"].upper(),
            "target": row["corrected_editing_efficiency"],
            "center_index": row["center_index"],
            "split": row["split"],
            "overlap_cluster_id": row["overlap_cluster_id"],
            "duplicate_sequence_group_id": row["duplicate_sequence_group_id"],
            "chrom": row["chrom"],
            "position": row["position"],
            "ref": row["ref"],
            "alt": row["alt"],
            "label_confidence": row["label_confidence"],
            "raw_edit_rate_difference": row["raw_edit_rate_difference"],
            "control_median": row["control_median"],
            "treated_median": row["treated_median"],
        }
        for row in rows
    ]


def validate_puf_target_requirement(token_mask_output, puf_target_seq):
    if token_mask_output is not None and not puf_target_seq:
        raise ValueError(
            "--token-mask-output requires the experimentally confirmed --puf-target-seq; "
            "the repository will not fabricate PUF metadata"
        )
    if puf_target_seq:
        sequence = puf_target_seq.upper()
        if any(base not in "ACGTU" for base in sequence):
            raise ValueError("PUF target sequence must contain only A/C/G/T/U")
        return sequence
    return None


def token_mask_rows(rows, puf_target_seq):
    output = []
    for row in rows:
        center = int(row["center_index"])
        labels = [None] * 101
        mask = [0] * 101
        labels[center] = float(row["corrected_editing_efficiency"])
        mask[center] = 1
        output.append(
            {
                "sequence": row["sequence"].upper(),
                "puf_target_seq": puf_target_seq,
                "label_values": json.dumps(labels, separators=(",", ":")),
                "label_mask": json.dumps(mask, separators=(",", ":")),
                "center_index": center,
                "split": row["split"],
                "chrom": row["chrom"],
                "position": row["position"],
                "ref": row["ref"],
                "alt": row["alt"],
            }
        )
    return output


def export_scalar(input_path, output_path, subset="high_confidence", token_mask_output=None, puf_target_seq=None):
    puf_target_seq = validate_puf_target_requirement(token_mask_output, puf_target_seq)
    rows = select_rows(read_tsv(input_path), subset)
    write_tsv(output_path, OUTPUT_FIELDS, scalar_rows(rows))
    if token_mask_output is not None:
        token_fields = (
            "sequence",
            "puf_target_seq",
            "label_values",
            "label_mask",
            "center_index",
            "split",
            "chrom",
            "position",
            "ref",
            "alt",
        )
        write_tsv(token_mask_output, token_fields, token_mask_rows(rows, puf_target_seq))
    return len(rows)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="CU5.17_lamar_splits.tsv.gz")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--subset",
        choices=("high_confidence", "all_eligible", "high_confidence_low_control"),
        default="high_confidence",
    )
    parser.add_argument("--token-mask-output", type=Path)
    parser.add_argument(
        "--puf-target-seq",
        help="Experimentally confirmed PUF target; mandatory only for optional token-mask export",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        count = export_scalar(
            args.input,
            args.output,
            subset=args.subset,
            token_mask_output=args.token_mask_output,
            puf_target_seq=args.puf_target_seq,
        )
    except Exception as error:
        raise SystemExit("ERROR: {}".format(error))
    print("exported_rows\t{}".format(count))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
