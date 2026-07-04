import argparse
import gzip
import os
import re
import sys
from pathlib import Path


def load_contig_order(fai_path):
    order = {}
    with open(fai_path, "r", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if not line.strip():
                continue
            chrom = line.split("\t", 1)[0]
            if chrom in order:
                raise ValueError("Duplicate contig in FAI: {}".format(chrom))
            order[chrom] = index
    if not order:
        raise ValueError("Reference FAI is empty: {}".format(fai_path))
    return order


def normalize_to_fai(chrom, order):
    if chrom in order:
        return chrom
    if "chr" + chrom in order:
        return "chr" + chrom
    if chrom.startswith("chr") and chrom[3:] in order:
        return chrom[3:]
    raise ValueError("Contig {!r} is not present in the reference FAI".format(chrom))


def expected_interval_count(temp_dir):
    path = temp_dir / "intervals.txt"
    if not path.is_file():
        raise ValueError("Missing intervals.txt: {}".format(path))
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = re.search(r"\d+", line)
            if match:
                return int(match.group(0))
    raise ValueError("Could not read expected interval count from {}".format(path))


def parse_chunk(path, order):
    name = path.name
    if not name.endswith(".gz"):
        return None

    # Remove only the final compression extension. Version suffixes in
    # contig identifiers, such as GL000009.2 and KI270750.1, must remain.
    stem = name[:-3]
    parts = stem.rsplit("#", 2)
    if len(parts) != 3:
        return None

    chrom, start_text, end_text = parts
    try:
        start = int(start_text)
        end = int(end_text)
    except ValueError:
        return None

    fai_chrom = normalize_to_fai(chrom, order)
    if path.stat().st_size == 0:
        raise ValueError("Empty REDItools2 interval file: {}".format(path))

    # Opening the stream catches truncated gzip files without reading the
    # entire interval into memory.
    try:
        with gzip.open(str(path), "rb") as handle:
            handle.read(1)
    except OSError as exc:
        raise ValueError("Invalid gzip interval file {}: {}".format(path, exc))

    return order[fai_chrom], start, end, str(path.resolve())


def rebuild(temp_dir, fai_path, output_path):
    order = load_contig_order(fai_path)
    expected = expected_interval_count(temp_dir)
    records = []
    seen = set()

    for path in temp_dir.glob("*.gz"):
        parsed = parse_chunk(path, order)
        if parsed is None:
            continue
        key = parsed[:3]
        if key in seen:
            raise ValueError("Duplicate REDItools2 interval: {}".format(key))
        seen.add(key)
        records.append(parsed)

    if len(records) != expected:
        raise ValueError(
            "Incomplete REDItools2 temporary output: expected {} interval files, found {}".format(
                expected, len(records)
            )
        )

    records.sort(key=lambda item: (item[0], item[1], item[2]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(output_path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for _contig_index, _start, _end, path in records:
            handle.write(path + "\n")
    os.replace(str(temporary), str(output_path))

    print("Expected intervals: {}".format(expected))
    print("Recovered intervals: {}".format(len(records)))
    print("File list: {}".format(output_path))
    if records:
        print("First chunk: {}".format(records[0][3]))
        print("Last chunk: {}".format(records[-1][3]))


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild a REDItools2 files.txt list while preserving contig version suffixes."
    )
    parser.add_argument("--temp-dir", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    temp_dir = Path(args.temp_dir).resolve()
    fai_path = Path(args.fai).resolve()
    output_path = Path(args.output).resolve() if args.output else temp_dir / "files.txt"

    try:
        rebuild(temp_dir, fai_path, output_path)
    except (OSError, ValueError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
