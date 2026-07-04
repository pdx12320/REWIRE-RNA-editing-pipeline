import argparse
import subprocess
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--cache", required=True)
    args = ap.parse_args()

    project = Path(args.project).resolve()
    vcf_path = project / "vcf" / "CU5.17_EGFP_GC.REDItools_union.vcf"
    output = project / "vep" / "CU5.17_EGFP_GC.vep.tsv"
    log_path = project / "logs" / "VEP.log"
    output.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "vep", "--cache", "--offline", "--dir_cache", str(Path(args.cache).resolve()),
        "--species", "homo_sapiens", "--assembly", "GRCh38", "--tab",
        "--fields", "Uploaded_variation,Location,Allele,STRAND",
        "--no_stats", "--input_file", str(vcf_path),
        "--output_file", str(output), "--force_overwrite",
    ]
    with open(log_path, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    if not output.exists() or output.stat().st_size == 0:
        raise SystemExit(f"VEP output is empty: {output}")
    print(output)


if __name__ == "__main__":
    main()
