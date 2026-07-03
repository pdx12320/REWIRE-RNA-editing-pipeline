import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

BASES = set("ACGT")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def chrom_key(chrom):
    x = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if x.isdigit():
        return (0, int(x))
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if x.upper() in order:
        return (0, order[x.upper()])
    return (1, x)


def parse_substitutions(raw, ref):
    if not raw or raw == "-":
        return []
    tokens = re.split(r"[\s,;]+", raw.strip())
    out = []
    for token in tokens:
        token = token.upper().replace(">", "")
        if len(token) == 2 and token[0] in BASES and token[1] in BASES:
            r, a = token
            if r == ref and a != r:
                out.append((r, a))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--reditools-dir", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bed", required=True)
    args = ap.parse_args()

    samples = []
    with open(args.manifest) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            samples.append(row["sample"])

    variants = defaultdict(set)
    for sample in samples:
        path = Path(args.reditools_dir) / f"{sample}.txt.gz"
        if not path.exists():
            raise SystemExit(f"Missing REDItools table: {path}")
        with open_text(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"Region", "Position", "Reference", "AllSubs"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise SystemExit(f"{path} missing REDItools columns: {sorted(missing)}")
            for row in reader:
                chrom = row["Region"]
                pos = int(row["Position"])
                ref = row["Reference"].upper()
                for r, alt in parse_substitutions(row["AllSubs"], ref):
                    variants[(chrom, pos, r, alt)].add(sample)

    ordered = sorted(variants, key=lambda x: (chrom_key(x[0]), x[1], x[2], x[3]))
    with open(args.vcf, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=REDItools2_union\n")
        out.write('##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Samples supporting this substitution">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in ordered:
            safe = re.sub(r"[^A-Za-z0-9_.-]", "_", chrom)
            vid = f"REDI_{safe}_{pos}_{ref}_{alt}"
            supp = ",".join(sorted(variants[(chrom, pos, ref, alt)]))
            out.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\tSAMPLES={supp}\n")

    positions = sorted({(c, p) for c, p, _, _ in ordered}, key=lambda x: (chrom_key(x[0]), x[1]))
    with open(args.bed, "w") as out:
        for chrom, pos in positions:
            out.write(f"{chrom}\t{pos-1}\t{pos}\n")

    print(f"variants={len(ordered)}")
    print(f"positions={len(positions)}")
    print(f"vcf={args.vcf}")
    print(f"bed={args.bed}")


if __name__ == "__main__":
    main()
