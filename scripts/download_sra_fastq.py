import argparse
import csv
import subprocess
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--manifest", default="config/samples.tsv")
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    project = Path(args.project).resolve()
    sra_dir = project / "sra"
    fastq_dir = project / "fastq"
    temp_dir = project / "tmp"
    for directory in [sra_dir, fastq_dir, temp_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    with open(args.manifest) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))

    for meta in samples:
        sample, srr = meta["sample"], meta["srr"]
        r1 = fastq_dir / f"{sample}_1.fastq.gz"
        r2 = fastq_dir / f"{sample}_2.fastq.gz"
        if r1.exists() and r2.exists():
            continue
        subprocess.run(["prefetch", srr, "--output-directory", str(sra_dir)], check=True)
        sra_file = sra_dir / srr / f"{srr}.sra"
        subprocess.run([
            "fasterq-dump", str(sra_file), "--split-files",
            "--threads", str(args.threads), "--temp", str(temp_dir),
            "--outdir", str(fastq_dir),
        ], check=True)
        for read in [1, 2]:
            source = fastq_dir / f"{srr}_{read}.fastq"
            target = fastq_dir / f"{sample}_{read}.fastq"
            source.rename(target)
            subprocess.run(["gzip", "-f", str(target)], check=True)


if __name__ == "__main__":
    main()
