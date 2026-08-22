#!/usr/bin/env python3
"""Audit 101-nt sequence uniqueness against the augmented reference with BWA."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
KEY_FIELDS = ("chrom", "position", "ref", "alt")


def open_text(path: Path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_rows(path: Path, sequence_column: str, eligible_only: bool = False):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input table has no header")
        missing = [field for field in KEY_FIELDS + (sequence_column,) if field not in reader.fieldnames]
        if missing:
            raise ValueError("Input table is missing columns: {}".format(", ".join(missing)))
        rows = []
        for row in reader:
            if eligible_only and str(row.get("training_eligible", "0")).strip().lower() not in {"1", "true", "yes", "y"}:
                continue
            sequence = row[sequence_column].upper().replace("U", "T")
            if not sequence or set(sequence) - set("ACGTN"):
                raise ValueError("Invalid DNA sequence for {}:{}".format(row["chrom"], row["position"]))
            row = dict(row)
            row[sequence_column] = sequence
            rows.append(row)
    if not rows:
        raise ValueError("Input table has no rows")
    return list(reader.fieldnames), rows


def full_length_alignment(cigar: str, query_length: int):
    operations = CIGAR_RE.findall(cigar)
    if not operations or "".join(count + op for count, op in operations) != cigar:
        return False
    query_consumed = sum(int(count) for count, op in operations if op in "MIS=X")
    clipped = any(op in "SH" for _, op in operations)
    return query_consumed == query_length and not clipped


def parse_nm(fields):
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            return int(field[5:])
    return None


def run_bwa(reference: Path, sequences, threads: int, bwa: str):
    stats = defaultdict(lambda: {"exact": 0, "nm2": 0, "mapped": 0, "best_nm": None})
    with tempfile.TemporaryDirectory(prefix="rewire_mappability_") as temporary:
        query_path = Path(temporary) / "queries.fa"
        with query_path.open("w", encoding="utf-8") as handle:
            for index, sequence in enumerate(sequences):
                handle.write(">q{}\n{}\n".format(index, sequence))
        command = [bwa, "mem", "-a", "-T", "0", "-t", str(threads), str(reference), str(query_path)]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        assert process.stdout is not None
        for line in process.stdout:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            query = fields[0]
            flag = int(fields[1])
            if flag & 0x4:
                continue
            index = int(query[1:])
            sequence = sequences[index]
            if not full_length_alignment(fields[5], len(sequence)):
                continue
            nm = parse_nm(fields)
            if nm is None:
                continue
            entry = stats[index]
            entry["mapped"] += 1
            entry["best_nm"] = nm if entry["best_nm"] is None else min(entry["best_nm"], nm)
            if nm == 0:
                entry["exact"] += 1
            if nm <= 2:
                entry["nm2"] += 1
        stderr = process.stderr.read() if process.stderr is not None else ""
        return_code = process.wait()
        if return_code:
            raise RuntimeError("BWA mappability audit failed: {}".format(stderr.strip()))
    return stats


def audit(input_path: Path, reference: Path, output: Path, sequence_column: str, threads: int, bwa: str, eligible_only: bool = False):
    fields, rows = read_rows(input_path, sequence_column, eligible_only=eligible_only)
    unique_sequences = sorted({row[sequence_column] for row in rows})
    index_by_sequence = {sequence: index for index, sequence in enumerate(unique_sequences)}
    stats = run_bwa(reference, unique_sequences, threads, bwa)
    output.parent.mkdir(parents=True, exist_ok=True)
    audit_fields = fields + [
        "exact_match_count",
        "nm_le_2_match_count",
        "full_length_mapped_count",
        "best_nm",
        "mappability_pass",
    ]
    passed = 0
    with open_text(output, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            entry = stats[index_by_sequence[row[sequence_column]]]
            is_pass = entry["exact"] == 1 and entry["nm2"] == 1
            passed += int(is_pass)
            writer.writerow({
                **row,
                "exact_match_count": entry["exact"],
                "nm_le_2_match_count": entry["nm2"],
                "full_length_mapped_count": entry["mapped"],
                "best_nm": entry["best_nm"] if entry["best_nm"] is not None else "NA",
                "mappability_pass": int(is_pass),
            })
    version_run = subprocess.run([bwa], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    version_lines = version_run.stdout.splitlines()
    version = version_lines[0] if version_lines else "unknown"
    report = {
        "status": "pass",
        "input": str(input_path.resolve()),
        "input_sha256": sha256(input_path),
        "reference": str(reference.resolve()),
        "reference_sha256": sha256(reference),
        "bwa_version_line": version,
        "rows": len(rows),
        "unique_sequences": len(unique_sequences),
        "mappability_pass_rows": passed,
        "mappability_fail_rows": len(rows) - passed,
        "pass_rule": "exact_match_count == 1 and nm_le_2_match_count == 1",
        "eligible_only": eligible_only,
    }
    report_path = Path(str(output) + ".audit.json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path, help="Augmented reference FASTA with BWA index")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sequence-column", default="sequence_context")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--bwa", default="bwa")
    parser.add_argument("--eligible-only", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        report = audit(
            args.input,
            args.reference,
            args.output,
            args.sequence_column,
            args.threads,
            args.bwa,
            eligible_only=args.eligible_only,
        )
    except Exception as error:
        raise SystemExit("ERROR: {}".format(error))
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
