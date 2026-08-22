#!/usr/bin/env python3
"""Append a verified EGFP-GC reporter contig to a genomic reference FASTA.

The paper identifies EGFP nucleotide 459 as the edited C and states that the
native U458-C459 context was changed to G458-C459.  It does not publish the
complete reporter FASTA.  This utility therefore requires a user-supplied,
sequence-verified reporter and refuses to guess it.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


DNA = set("ACGTN")


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")


def read_fasta(path: Path):
    name = None
    chunks = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                if not name:
                    raise ValueError("FASTA contains an empty record name")
                chunks = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence appears before the first header")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def wrap(sequence: str, width: int = 60):
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_reporter(sequence: str, target_position: int, expected_context: str):
    if not sequence:
        raise ValueError("Reporter FASTA sequence is empty")
    invalid = sorted(set(sequence) - DNA)
    if invalid:
        raise ValueError("Reporter FASTA contains non-DNA symbols: {}".format("".join(invalid)))
    expected_context = expected_context.upper()
    if len(expected_context) != 2 or set(expected_context) - DNA:
        raise ValueError("Expected context must contain exactly two DNA bases")
    if target_position < 2 or target_position > len(sequence):
        raise ValueError("Target position is outside the reporter sequence")
    observed = sequence[target_position - 2 : target_position]
    if observed != expected_context:
        raise ValueError(
            "Reporter validation failed at 1-based positions {}-{}: expected {}, observed {}".format(
                target_position - 1, target_position, expected_context, observed
            )
        )
    if sequence[target_position - 1] != "C":
        raise AssertionError("Validated target base is not C")
    return observed


def build_reference(
    genome_fasta: Path,
    reporter_fasta: Path,
    output_fasta: Path,
    contig_name: str,
    target_position: int,
    expected_context: str,
    output_gtf: Path | None = None,
    genome_gtf: Path | None = None,
):
    reporter_records = list(read_fasta(reporter_fasta))
    if len(reporter_records) != 1:
        raise ValueError("Reporter FASTA must contain exactly one record")
    source_name, reporter_sequence = reporter_records[0]
    validate_reporter(reporter_sequence, target_position, expected_context)

    genome_names = []
    for name, _ in read_fasta(genome_fasta):
        genome_names.append(name)
    if contig_name in genome_names:
        raise ValueError("Reporter contig already exists in genome FASTA: {}".format(contig_name))

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as out:
        ended_with_newline = True
        with open_text(genome_fasta) as source:
            for line in source:
                out.write(line)
                ended_with_newline = line.endswith("\n")
        if not ended_with_newline:
            out.write("\n")
        out.write(">{} source={} target=C{} context={}\n".format(
            contig_name, source_name, target_position, expected_context.upper()
        ))
        for line in wrap(reporter_sequence):
            out.write(line + "\n")

    if output_gtf is not None:
        if genome_gtf is None:
            raise ValueError("--genome-gtf is required when --output-gtf is used")
        output_gtf.parent.mkdir(parents=True, exist_ok=True)
        with output_gtf.open("w", encoding="utf-8") as out:
            with open_text(genome_gtf) as source:
                for line in source:
                    out.write(line)
            attributes = 'gene_id "EGFP_GC_reporter"; gene_name "EGFP_GC_reporter"; transcript_id "EGFP_GC_reporter";'
            out.write("{}\tREWIRE\texon\t1\t{}\t.\t+\t.\t{}\n".format(
                contig_name, len(reporter_sequence), attributes
            ))

    report = {
        "genome_fasta": str(genome_fasta.resolve()),
        "genome_sha256": sha256(genome_fasta),
        "reporter_fasta": str(reporter_fasta.resolve()),
        "reporter_sha256": sha256(reporter_fasta),
        "source_reporter_record": source_name,
        "output_fasta": str(output_fasta.resolve()),
        "output_sha256": sha256(output_fasta),
        "reporter_contig": contig_name,
        "reporter_length": len(reporter_sequence),
        "target_position_1based": target_position,
        "target_base": reporter_sequence[target_position - 1],
        "target_context": reporter_sequence[target_position - 2 : target_position],
        "status": "pass",
    }
    audit_path = output_fasta.with_suffix(output_fasta.suffix + ".audit.json")
    audit_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-fasta", required=True, type=Path)
    parser.add_argument("--reporter-fasta", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--contig-name", default="EGFP_GC_reporter")
    parser.add_argument("--target-position", type=int, default=459)
    parser.add_argument("--expected-context", default="GC")
    parser.add_argument("--genome-gtf", type=Path)
    parser.add_argument("--output-gtf", type=Path)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.output_gtf is not None and args.genome_gtf is None:
        raise SystemExit("ERROR: --genome-gtf is required with --output-gtf")
    try:
        report = build_reference(
            args.genome_fasta,
            args.reporter_fasta,
            args.output_fasta,
            args.contig_name,
            args.target_position,
            args.expected_context,
            args.output_gtf,
            args.genome_gtf,
        )
    except Exception as error:
        raise SystemExit("ERROR: {}".format(error))
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
