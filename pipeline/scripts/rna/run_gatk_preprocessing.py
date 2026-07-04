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
    ap.add_argument("--reference", required=True)
    ap.add_argument("--manifest", default="config/samples.tsv")
    ap.add_argument("--java-options", default="-Xmx16g")
    args = ap.parse_args()

    project = Path(args.project).resolve()
    star_dir = project / "bam" / "star"
    rg_dir = project / "bam" / "readgroups"
    mark_dir = project / "bam" / "markduplicates"
    split_dir = project / "bam" / "splitncigarreads"
    metrics_dir = project / "metrics"
    log_dir = project / "logs"
    for directory in [rg_dir, mark_dir, split_dir, metrics_dir, log_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    with open(args.manifest) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))

    for meta in samples:
        sample = meta["sample"]
        star = star_dir / f"{sample}.Aligned.sortedByCoord.out.bam"
        rg = rg_dir / f"{sample}.readgroups.bam"
        mark = mark_dir / f"{sample}.markduplicates.bam"
        split = split_dir / f"{sample}.splitncigarreads.bam"
        metrics = metrics_dir / f"{sample}.markduplicates.txt"
        if not star.exists():
            raise SystemExit(f"Missing STAR BAM: {star}")

        header = subprocess.check_output(["samtools", "view", "-H", str(star)], text=True)
        source = star
        if not any(line.startswith("@RG") for line in header.splitlines()):
            if not rg.exists() or rg.stat().st_size == 0:
                cmd = [
                    "gatk", "--java-options", args.java_options,
                    "AddOrReplaceReadGroups", "-I", str(star), "-O", str(rg),
                    "--RGID", sample, "--RGLB", sample, "--RGPL", "ILLUMINA",
                    "--RGPU", sample, "--RGSM", sample,
                ]
                run(cmd, log_dir / f"{sample}.readgroups.log")
            source = rg

        if not mark.exists() or mark.stat().st_size == 0:
            cmd = [
                "gatk", "--java-options", args.java_options,
                "MarkDuplicates", "-I", str(source), "-O", str(mark),
                "-M", str(metrics), "--CREATE_INDEX", "true",
            ]
            run(cmd, log_dir / f"{sample}.markduplicates.log")

        if not split.exists() or split.stat().st_size == 0:
            cmd = [
                "gatk", "--java-options", args.java_options,
                "SplitNCigarReads", "-R", str(Path(args.reference).resolve()),
                "-I", str(mark), "-O", str(split),
            ]
            run(cmd, log_dir / f"{sample}.splitncigarreads.log")
            subprocess.run(["samtools", "index", str(split)], check=True)
        subprocess.run(["samtools", "quickcheck", "-v", str(split)], check=True)


if __name__ == "__main__":
    main()
