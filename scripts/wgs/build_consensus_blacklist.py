#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fai(path):
    order = {}
    lengths = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            p = line.split("\t")
            order[p[0]] = i
            lengths[p[0]] = int(p[1])
    return order, lengths


def main():
    ap = argparse.ArgumentParser(description="Build exact-allele WGS union and N-of-M consensus VCFs.")
    ap.add_argument("--vcf", action="append", required=True, help="Filtered single-sample VCF; repeat for each run")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--min-support", type=int, default=2)
    ap.add_argument("--union-output", required=True)
    ap.add_argument("--consensus-output", required=True)
    args = ap.parse_args()

    order, lengths = read_fai(args.fai)
    counts = defaultdict(set)
    records = {}
    for idx, path in enumerate(args.vcf):
        seen = set()
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 8:
                    continue
                chrom, pos, vid, ref, alt, qual, filt, info = p[:8]
                key = (chrom, int(pos), ref, alt)
                if key in seen:
                    continue
                seen.add(key)
                counts[key].add(idx)
                records.setdefault(key, (qual, filt, info))

    def write(path, min_support):
        with open(path, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("##source=REWIRE_public_HEK293T_WGS_blacklist\n")
            out.write('##INFO=<ID=RUN_SUPPORT,Number=1,Type=Integer,Description="Number of input WGS runs supporting the exact CHROM:POS:REF:ALT allele">\n')
            for chrom in order:
                out.write("##contig=<ID={},length={}>\n".format(chrom, lengths[chrom]))
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            keys = [k for k, v in counts.items() if len(v) >= min_support]
            keys.sort(key=lambda x: (order.get(x[0], 10**9), x[1], x[2], x[3]))
            for chrom, pos, ref, alt in keys:
                qual, filt, info = records[(chrom, pos, ref, alt)]
                support = len(counts[(chrom, pos, ref, alt)])
                new_info = ("" if info in {"", ".", "-"} else info + ";") + "RUN_SUPPORT={}".format(support)
                out.write("{}\t{}\t.\t{}\t{}\t{}\t{}\t{}\n".format(chrom, pos, ref, alt, qual, filt, new_info))
        return len(keys)

    nu = write(args.union_output, 1)
    nc = write(args.consensus_output, args.min_support)
    print("Union variants:", nu)
    print("Consensus variants:", nc)


if __name__ == "__main__":
    main()
