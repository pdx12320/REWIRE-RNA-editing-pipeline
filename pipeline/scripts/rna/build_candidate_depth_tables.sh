set -Eeuo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: bash $0 PROJECT UNION_BED MANIFEST" >&2
  exit 2
fi

PROJECT="$(readlink -m "$1")"
UNION_BED="$(readlink -m "$2")"
MANIFEST="$(readlink -m "$3")"
BAM_DIR="$PROJECT/bam/splitncigarreads"
OUT_DIR="$PROJECT/candidate_depth"
mkdir -p "$OUT_DIR"

while IFS=$'\t' read -r sample group replicate srr; do
  [[ "$sample" == "sample" ]] && continue
  bam="$BAM_DIR/${sample}.splitncigarreads.bam"
  out="$OUT_DIR/${sample}.candidate_depth.tsv.gz"
  [[ -f "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }
  samtools depth -a -q 30 -Q 20 -b "$UNION_BED" "$bam" | gzip -c > "$out"
  [[ -s "$out" ]] || { echo "Empty depth table: $out" >&2; exit 1; }
  echo "Candidate depth complete: $sample" >&2
done < "$MANIFEST"
