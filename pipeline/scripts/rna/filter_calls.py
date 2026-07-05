import csv
from collections import defaultdict
from pathlib import Path

from filter_utils import BASE_INDEX, norm_chrom, open_text, parse_counts, parse_substitutions, variant_id

CALL_FIELDS = [
    "sample", "group", "replicate", "srr", "chrom", "position", "ref", "alt",
    "vep_strand", "coverage_q30", "mean_q", "ref_count", "alt_count",
    "edit_rate", "all_substitutions", "genomic_catalogue_overlap"
]


def collect_calls(manifest, reditools_dir, output_dir, vep, catalogue_variants):
    all_calls = defaultdict(dict)
    sample_dir = Path(output_dir) / "sample_calls"
    sample_dir.mkdir(parents=True, exist_ok=True)

    for meta in manifest:
        sample = meta["sample"]
        src = Path(reditools_dir) / (sample + ".txt.gz")
        dst = sample_dir / (sample + ".c_to_u.tsv.gz")
        with open_text(src) as fh, open_text(dst, "wt") as out:
            reader = csv.DictReader(fh, delimiter="\t")
            writer = csv.DictWriter(out, fieldnames=CALL_FIELDS, delimiter="\t")
            writer.writeheader()
            for row in reader:
                chrom, pos = row["Region"], int(row["Position"])
                ref = row["Reference"].upper()
                counts = parse_counts(row["BaseCount[A,C,G,T]"])
                for r, alt in parse_substitutions(row["AllSubs"], ref):
                    vid = variant_id(chrom, pos, r, alt)
                    strand = vep.get((vid, alt), vep.get((vid, "")))
                    is_c_to_u = (
                        strand == 1 and r == "C" and alt == "T"
                    ) or (
                        strand == -1 and r == "G" and alt == "A"
                    )
                    if not is_c_to_u:
                        continue
                    ref_count = counts[BASE_INDEX[r]]
                    alt_count = counts[BASE_INDEX[alt]]
                    denom = ref_count + alt_count
                    rec = {
                        "sample": sample,
                        "group": meta["group"],
                        "replicate": meta["replicate"],
                        "srr": meta["srr"],
                        "chrom": chrom,
                        "position": pos,
                        "ref": r,
                        "alt": alt,
                        "vep_strand": strand,
                        "coverage_q30": int(float(row["Coverage-q30"])),
                        "mean_q": row.get("MeanQ", ""),
                        "ref_count": ref_count,
                        "alt_count": alt_count,
                        "edit_rate": "{:.8g}".format(alt_count / float(denom) if denom else 0.0),
                        "all_substitutions": row["AllSubs"],
                        "genomic_catalogue_overlap": int(
                            (norm_chrom(chrom), pos, r, alt) in catalogue_variants
                        ),
                    }
                    writer.writerow(rec)
                    all_calls[(chrom, pos, r, alt)][sample] = rec
    return all_calls
