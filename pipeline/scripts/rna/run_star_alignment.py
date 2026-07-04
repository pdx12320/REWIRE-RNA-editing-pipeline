import argparse
import csv
import subprocess
from pathlib import Path


def run(cmd, log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--star-index", required=True)
    ap.add_argument("--manifest", default="config/samples.tsv")
    ap.add_argument("--threads", type=int, default=50)
    args = ap.parse_args()

    project = Path(args.project).resolve()
    fastq_dir = project / "fastq"
    bam_dir = project / "bam" / "star"
    log_dir = project / "logs"
    bam_dir.mkdir(parents=True, exist_ok=True)

    with open(args.manifest) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))

    for meta in samples:
        sample = meta["sample"]
        r1 = fastq_dir / f"{sample}_1.fastq.gz"
        r2 = fastq_dir / f"{sample}_2.fastq.gz"
        output = bam_dir / f"{sample}.Aligned.sortedByCoord.out.bam"
        if output.exists() and output.stat().st_size > 0:
            print("skip", sample)
            continue
        if not r1.exists() or not r2.exists():
            raise SystemExit(f"Missing FASTQ files for {sample}")
        prefix = str(bam_dir / f"{sample}.")
        cmd = [
            "STAR", "--genomeDir", str(Path(args.star_index).resolve()),
            "--runThreadN", str(args.threads),
            "--readFilesIn", str(r1), str(r2),
            "--readFilesCommand", "gunzip", "-c",
            "--twopassMode", "Basic",
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMattrRGline", f"ID:{sample}", f"SM:{sample}",
            "PL:ILLUMINA", f"LB:{sample}", f"PU:{sample}",
            "--outFileNamePrefix", prefix,
        ]
        run(cmd, log_dir / f"{sample}.STAR.log")
        subprocess.run(["samtools", "index", "-@", str(args.threads), str(output)], check=True)
        subprocess.run(["samtools", "quickcheck", "-v", str(output)], check=True)


if __name__ == "__main__":
    main()
