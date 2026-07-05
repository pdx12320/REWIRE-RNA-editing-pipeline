import argparse
import csv
from pathlib import Path

from filter_calls import collect_calls
from filter_utils import load_depth, norm_chrom, open_text, parse_variant_catalogue_vcf, parse_vep


def main():
    ap = argparse.ArgumentParser(
        description="Integrate strand-aware REDItools2 calls, all-sample depth and an optional genomic variant catalogue."
    )
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--reditools-dir", required=True)
    ap.add_argument("--vep", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--depth-dir", required=True)
    ap.add_argument(
        "--variant-catalogue-vcf",
        default="",
        help="Normalized exact-allele VCF used as an external genomic-variant catalogue.",
    )
    ap.add_argument(
        "--wgs-vcf",
        default="",
        help="Deprecated compatibility alias for --variant-catalogue-vcf.",
    )
    ap.add_argument("--min-treated-reps", type=int, default=3)
    ap.add_argument("--max-control-called-reps", type=int, default=0)
    ap.add_argument("--min-depth-all-reps", type=int, default=20)
    args = ap.parse_args()

    if args.variant_catalogue_vcf and args.wgs_vcf:
        ap.error("use only one of --variant-catalogue-vcf and --wgs-vcf")
    catalogue_path = args.variant_catalogue_vcf or args.wgs_vcf
    if args.wgs_vcf:
        print("WARNING: --wgs-vcf is deprecated; use --variant-catalogue-vcf")

    with open(args.manifest) as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    samples = [m["sample"] for m in manifest]
    treated = [m["sample"] for m in manifest if m["group"] == "treated"]
    controls = [m["sample"] for m in manifest if m["group"] == "control"]

    depths = {}
    for sample in samples:
        path = Path(args.depth_dir) / (sample + ".candidate_depth.tsv.gz")
        depths[sample] = load_depth(path)

    calls = collect_calls(
        manifest,
        args.reditools_dir,
        args.output_dir,
        parse_vep(args.vep),
        parse_variant_catalogue_vcf(catalogue_path),
    )

    fields = [
        "chrom", "position", "ref", "alt", "vep_strand",
        "treated_called_reps", "control_called_reps",
        "all_replicates_depth_pass", "minimum_replicate_depth",
        "genomic_catalogue_overlap"
    ]
    for sample in samples:
        fields += [
            sample + "_called",
            sample + "_candidate_depth_q30",
            sample + "_reditools_coverage_q30",
            sample + "_alt_count",
            sample + "_edit_rate",
        ]

    rows = []
    for (chrom, pos, ref, alt), per_sample in sorted(calls.items()):
        d = {s: depths[s].get((norm_chrom(chrom), pos), 0) for s in samples}
        row = {
            "chrom": chrom,
            "position": pos,
            "ref": ref,
            "alt": alt,
            "vep_strand": next(iter(per_sample.values()))["vep_strand"],
            "treated_called_reps": sum(s in per_sample for s in treated),
            "control_called_reps": sum(s in per_sample for s in controls),
            "all_replicates_depth_pass": int(
                all(v >= args.min_depth_all_reps for v in d.values())
            ),
            "minimum_replicate_depth": min(d.values()) if d else 0,
            "genomic_catalogue_overlap": int(
                any(bool(x["genomic_catalogue_overlap"]) for x in per_sample.values())
            ),
        }
        for sample in samples:
            rec = per_sample.get(sample)
            row[sample + "_called"] = int(rec is not None)
            row[sample + "_candidate_depth_q30"] = d[sample]
            row[sample + "_reditools_coverage_q30"] = rec["coverage_q30"] if rec else "NA"
            row[sample + "_alt_count"] = rec["alt_count"] if rec else "NA"
            row[sample + "_edit_rate"] = rec["edit_rate"] if rec else "NA"
        rows.append(row)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    consensus = [
        r for r in rows
        if int(r["treated_called_reps"]) >= args.min_treated_reps
        and int(r["all_replicates_depth_pass"]) == 1
    ]
    specific = [
        r for r in consensus
        if int(r["control_called_reps"]) <= args.max_control_called_reps
        and int(r["genomic_catalogue_overlap"]) == 0
    ]
    datasets = {
        "CU5.17_EGFP_GC.site_matrix.tsv.gz": rows,
        "CU5.17_EGFP_GC.treated_consensus.tsv.gz": consensus,
        "CU5.17_EGFP_GC.treatment_specific.tsv.gz": specific,
    }

    for name, data in datasets.items():
        with open_text(outdir / name, "wt") as out:
            writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(data)
        print(name, len(data))


if __name__ == "__main__":
    main()
