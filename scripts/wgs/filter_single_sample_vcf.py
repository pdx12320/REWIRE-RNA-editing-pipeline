#!/usr/bin/env python3
import argparse
import gzip


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser(description="Filter a normalized single-sample biallelic SNP VCF.")
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--min-alt", type=int, default=3)
    ap.add_argument("--min-vaf", type=float, default=0.05)
    ap.add_argument("--min-qual", type=float, default=20.0)
    args = ap.parse_args()

    kept = 0
    total = 0
    with open_text(args.input) as inp, open(args.output, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info, fmt, sample = fields[:10]
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            try:
                q = float(qual) if qual != "." else 0.0
            except ValueError:
                q = 0.0
            if q < args.min_qual or filt not in {".", "PASS"}:
                continue
            keys = fmt.split(":")
            vals = sample.split(":")
            data = dict(zip(keys, vals))
            try:
                dp = int(data.get("DP", "0") or 0)
            except ValueError:
                dp = 0
            ad_raw = data.get("AD", "")
            try:
                ad = [int(x) for x in ad_raw.split(",")]
            except Exception:
                ad = []
            if len(ad) < 2:
                continue
            ref_count, alt_count = ad[0], ad[1]
            denom = ref_count + alt_count
            vaf = alt_count / float(denom) if denom else 0.0
            if dp < args.min_dp or alt_count < args.min_alt or vaf < args.min_vaf:
                continue
            fields[7] = ("" if info == "." else info + ";") + "WGS_DP={};WGS_ALT={};WGS_VAF={:.6g}".format(dp, alt_count, vaf)
            out.write("\t".join(fields) + "\n")
            kept += 1
    print("Input variants:", total)
    print("Retained SNPs:", kept)


if __name__ == "__main__":
    main()
