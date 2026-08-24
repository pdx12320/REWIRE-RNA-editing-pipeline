#!/usr/bin/env python3
"""Audit that every sample uses the same analysis-ready BAM processing route."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


REQUIRED_PROGRAM_TOKENS = ("STAR", "MarkDuplicates", "SplitNCigarReads")


def read_manifest(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"sample", "group", "replicate"}.issubset(rows[0]):
        raise ValueError("Manifest must contain sample, group and replicate columns")
    names = [row["sample"] for row in rows]
    if len(names) != len(set(names)):
        raise ValueError("Manifest contains duplicate sample names")
    return rows


def dictionary_digest(alignment):
    text = "\n".join("{}\t{}".format(name, length) for name, length in zip(alignment.references, alignment.lengths))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def bam_path(bam_dir: Path, sample: str):
    path = bam_dir / "{}.splitncigarreads.bam".format(sample)
    if not path.is_file():
        raise FileNotFoundError("Missing analysis-ready BAM: {}".format(path))
    return path


def audit_one(path: Path, sample_reads: int, expected_contig: str | None):
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required for BAM auditing") from exc

    subprocess.run(["samtools", "quickcheck", "-v", str(path)], check=True)
    with pysam.AlignmentFile(path, "rb") as alignment:
        header = alignment.header.to_dict()
        program_text = json.dumps(header.get("PG", []), sort_keys=True)
        program_presence = {
            token: token.lower() in program_text.lower() for token in REQUIRED_PROGRAM_TOKENS
        }
        if not all(program_presence.values()):
            missing = [token for token, present in program_presence.items() if not present]
            raise ValueError("{} lacks required @PG provenance: {}".format(path, ", ".join(missing)))
        if expected_contig and expected_contig not in alignment.references:
            raise ValueError("{} lacks expected reporter contig {}".format(path, expected_contig))

        inspected = 0
        missing_nh = 0
        nh1 = 0
        nh_gt1 = 0
        for read in alignment.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            inspected += 1
            if not read.has_tag("NH"):
                missing_nh += 1
            elif read.get_tag("NH") == 1:
                nh1 += 1
            else:
                nh_gt1 += 1
            if inspected >= sample_reads:
                break
        if inspected == 0:
            raise ValueError("{} contains no mapped primary reads to audit".format(path))
        if missing_nh:
            raise ValueError("{} has {} / {} sampled primary reads without NH".format(path, missing_nh, inspected))
        return {
            "bam": str(path.resolve()),
            "reference_dictionary_sha256": dictionary_digest(alignment),
            "mapped_primary_reads_sampled": inspected,
            "missing_nh": missing_nh,
            "nh_eq_1": nh1,
            "nh_gt_1": nh_gt1,
            "nh_eq_1_fraction": nh1 / inspected,
            "required_programs": program_presence,
            "reporter_contig_present": expected_contig in alignment.references if expected_contig else None,
            "status": "pass",
        }


def audit_all(manifest: Path, bam_dir: Path, output: Path, expected_samples: int, sample_reads: int, expected_contig: str | None):
    rows = read_manifest(manifest)
    if expected_samples and len(rows) != expected_samples:
        raise ValueError("Expected {} samples, found {}".format(expected_samples, len(rows)))
    audits = []
    for row in rows:
        audit = audit_one(bam_path(bam_dir, row["sample"]), sample_reads, expected_contig)
        audit.update({"sample": row["sample"], "group": row["group"], "replicate": row["replicate"]})
        audits.append(audit)
    dictionaries = {audit["reference_dictionary_sha256"] for audit in audits}
    if len(dictionaries) != 1:
        raise ValueError("The six BAMs do not share one reference sequence dictionary")
    result = {
        "status": "pass",
        "sample_count": len(audits),
        "bam_stage": "splitncigarreads",
        "required_programs": list(REQUIRED_PROGRAM_TOKENS),
        "expected_reporter_contig": expected_contig,
        "reference_dictionary_sha256": next(iter(dictionaries)),
        "samples": audits,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--bam-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-samples", type=int, default=6)
    parser.add_argument("--sample-reads", type=int, default=100000)
    parser.add_argument("--expected-contig", default="EGFP_GC_reporter")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        result = audit_all(
            args.manifest,
            args.bam_dir,
            args.output,
            args.expected_samples,
            args.sample_reads,
            args.expected_contig or None,
        )
    except Exception as error:
        raise SystemExit("ERROR: {}".format(error))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
