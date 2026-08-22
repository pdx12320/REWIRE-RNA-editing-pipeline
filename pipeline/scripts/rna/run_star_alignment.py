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
    ap.add_argument(
        "--expected-contig",
        default=None,
        help="Abort unless this contig is present in every output BAM (use EGFP_GC_reporter for the strict route)",
    )
    ap.add_argument("--expected-samples", type=int, default=6)
    ap.add_argument("--force", action="store_true", help="Replace only this script's STAR BAM and index outputs")
    args = ap.parse_args()

    project = Path(args.project).resolve()
    fastq_dir = project / "fastq"
    bam_dir = project / "bam" / "star"
    log_dir = project / "logs"
    bam_dir.mkdir(parents=True, exist_ok=True)

    with open(args.manifest) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))
    if args.expected_samples and len(samples) != args.expected_samples:
        raise SystemExit(f"Expected {args.expected_samples} samples, found {len(samples)}")

    for meta in samples:
        sample = meta["sample"]
        r1 = fastq_dir / f"{sample}_1.fastq.gz"
        r2 = fastq_dir / f"{sample}_2.fastq.gz"
        output = bam_dir / f"{sample}.Aligned.sortedByCoord.out.bam"
        if args.force:
            for generated in (output, Path(str(output) + ".bai")):
                if generated.is_file():
                    generated.unlink()
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
            "--outSAMattributes", "NH", "HI", "AS", "nM", "MD",
            "--outSAMmapqUnique", "60",
            "--outSAMattrRGline", f"ID:{sample}", f"SM:{sample}",
            "PL:ILLUMINA", f"LB:{sample}", f"PU:{sample}",
            "--outFileNamePrefix", prefix,
        ]
        run(cmd, log_dir / f"{sample}.STAR.log")
        subprocess.run(["samtools", "index", "-@", str(args.threads), str(output)], check=True)
        subprocess.run(["samtools", "quickcheck", "-v", str(output)], check=True)
        if args.expected_contig:
            header = subprocess.check_output(["samtools", "view", "-H", str(output)], text=True)
            sq_names = {
                token[3:]
                for line in header.splitlines()
                if line.startswith("@SQ")
                for token in line.split("\t")
                if token.startswith("SN:")
            }
            if args.expected_contig not in sq_names:
                raise SystemExit(
                    f"Expected contig {args.expected_contig!r} is absent from {output}; rebuild STAR index with the augmented FASTA"
                )


if __name__ == "__main__":
    main()
