#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
cat <<'USAGE'
Convert the HEK293 Genome Project 293T_CG VCF from NCBI build 36 / hg18
into a GRCh38 exact-allele catalogue.

Usage:
  bash process_293T_CG_to_GRCh38.sh \
    --input /path/293T_CG.vcf[.gz] \
    --reference /path/GRCh38.fa \
    --chain /path/hg18ToHg38.over.chain.gz \
    --outdir /path/293T_CG_GRCh38 \
    [--threads 16]

Outputs include the full normalized catalogue and a C>T/G>A subset for
C-to-U RNA-editing evidence integration.
USAGE
}

INPUT=""
REF=""
CHAIN=""
OUTDIR=""
THREADS=16

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT="$2"; shift 2 ;;
    --reference) REF="$2"; shift 2 ;;
    --chain) CHAIN="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -s "$INPUT" ]] || { echo "ERROR: input VCF not found: $INPUT" >&2; exit 1; }
[[ -s "$REF" ]] || { echo "ERROR: GRCh38 FASTA not found: $REF" >&2; exit 1; }
[[ -s "$CHAIN" ]] || { echo "ERROR: hg18ToHg38 chain not found: $CHAIN" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 1; }

for cmd in bcftools bgzip tabix samtools python3 sha256sum; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing command: $cmd" >&2; exit 1; }
done
if command -v CrossMap.py >/dev/null 2>&1; then
  CROSSMAP="CrossMap.py"
elif command -v CrossMap >/dev/null 2>&1; then
  CROSSMAP="CrossMap"
else
  echo "ERROR: CrossMap is not installed" >&2
  exit 1
fi

gzip -t "$CHAIN"
mkdir -p "$OUTDIR"/{logs,qc,tmp}
LOG="$OUTDIR/logs/process.log"
exec > >(tee -a "$LOG") 2>&1
trap 'echo "[ERROR] Failed at line $LINENO. See: $LOG" >&2' ERR

echo "[$(date '+%F %T')] Starting 293T catalogue conversion"
echo "INPUT   = $INPUT"
echo "REF     = $REF"
echo "CHAIN   = $CHAIN"
echo "OUTDIR  = $OUTDIR"
echo "THREADS = $THREADS"

[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

# Normalize the input by content rather than extension. This handles a gzip
# stream saved as .vcf and files with harmless bytes before the gzip header.
NORMALIZED_INPUT="$OUTDIR/tmp/293T_CG.input.normalized.vcf"
python3 - "$INPUT" "$NORMALIZED_INPUT" <<'PY'
import gzip
import shutil
import sys

src, dst = sys.argv[1:3]
with open(src, "rb") as fh:
    prefix = fh.read(65536)

stripped = prefix.lstrip()
if stripped.startswith(b"##fileformat=VCF"):
    with open(src, "rb") as inp, open(dst, "wb") as out:
        shutil.copyfileobj(inp, out)
else:
    offset = prefix.find(b"\x1f\x8b\x08")
    if offset < 0:
        if b"<html" in prefix.lower() or b"<!doctype html" in prefix.lower():
            raise SystemExit("ERROR: input is an HTML page, not a VCF")
        raise SystemExit("ERROR: input is neither a recognizable text VCF nor gzip VCF")
    with open(src, "rb") as raw:
        raw.seek(offset)
        with gzip.GzipFile(fileobj=raw, mode="rb") as inp, open(dst, "wb") as out:
            shutil.copyfileobj(inp, out)

with open(dst, "rb") as fh:
    first = fh.readline().strip()
if not first.startswith(b"##fileformat=VCF"):
    raise SystemExit("ERROR: normalized input does not start with ##fileformat=VCF")
PY

HEADER="$OUTDIR/qc/input_header.txt"
grep '^#' "$NORMALIZED_INPUT" > "$HEADER"
if ! grep -Eqi 'build36|hg18|NCBI[ _-]?36' "$HEADER"; then
  echo "WARNING: input header does not explicitly identify build36/hg18" >&2
fi
grep -E '^##(fileformat|fileDate|source|reference)=' "$HEADER" || true

SOURCE_SNV="$OUTDIR/293T_CG.hg18.PASS.biallelic.SNV.vcf.gz"
echo "[$(date '+%F %T')] Selecting source PASS biallelic SNPs"
bcftools view \
  --threads "$THREADS" \
  -f PASS -v snps -m2 -M2 -i 'GT="alt"' \
  "$NORMALIZED_INPUT" \
  -Oz -o "$SOURCE_SNV"
bcftools index --threads "$THREADS" -f -t "$SOURCE_SNV"
SOURCE_N=$(bcftools index -n "$SOURCE_SNV")
echo "Source PASS biallelic SNPs: $SOURCE_N"
[[ "$SOURCE_N" -ge 1000000 ]] || echo "WARNING: fewer than one million source SNPs retained" >&2

LIFTED="$OUTDIR/tmp/293T_CG.GRCh38.CrossMap.vcf"
rm -f "$LIFTED" "${LIFTED}.unmap"
echo "[$(date '+%F %T')] Lifting hg18 coordinates to GRCh38"
"$CROSSMAP" vcf "$CHAIN" "$SOURCE_SNV" "$REF" "$LIFTED"
[[ -s "$LIFTED" ]] || { echo "ERROR: CrossMap produced no VCF" >&2; exit 1; }

UNMAPPED=0
if [[ -s "${LIFTED}.unmap" ]]; then
  UNMAPPED=$(grep -vc '^#' "${LIFTED}.unmap" || true)
fi
echo "CrossMap unmapped records: $UNMAPPED"

FINAL="$OUTDIR/293T_CG.GRCh38.PASS.biallelic.SNV.vcf.gz"
SORT_TMP="$OUTDIR/tmp/bcftools_sort"
rm -rf "$SORT_TMP"
mkdir -p "$SORT_TMP"
rm -f "$FINAL" "${FINAL}.tbi" "${FINAL}.csi"

echo "[$(date '+%F %T')] Validating REF, normalizing and sorting"
bcftools norm \
  --threads "$THREADS" \
  -f "$REF" -c x -m -any \
  "$LIFTED" -Ou | \
bcftools view \
  --threads "$THREADS" \
  -f PASS -v snps -m2 -M2 -i 'GT="alt"' \
  -Ou | \
bcftools sort \
  -T "$SORT_TMP" \
  -Oz -o "$FINAL"
bcftools index --threads "$THREADS" -f -t "$FINAL"

CTGA="$OUTDIR/293T_CG.GRCh38.CtoU_relevant.SNV.vcf.gz"
bcftools view \
  --threads "$THREADS" \
  -i '(REF="C" && ALT="T") || (REF="G" && ALT="A")' \
  "$FINAL" -Oz -o "$CTGA"
bcftools index --threads "$THREADS" -f -t "$CTGA"

FINAL_N=$(bcftools index -n "$FINAL")
CTGA_N=$(bcftools index -n "$CTGA")
MAPPED_N=$((SOURCE_N - UNMAPPED))
MISMATCH_REMOVED=$((MAPPED_N - FINAL_N))

bcftools stats "$FINAL" > "$OUTDIR/qc/293T_CG.GRCh38.stats.txt"
bcftools query -f '%CHROM\n' "$FINAL" | sort -V | uniq -c > "$OUTDIR/qc/variants_per_contig.txt"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$FINAL" | head -n 20 > "$OUTDIR/qc/first_20_variants.tsv"
sha256sum "$FINAL" "${FINAL}.tbi" "$CTGA" "${CTGA}.tbi" > "$OUTDIR/qc/checksums.sha256"

cat > "$OUTDIR/qc/summary.txt" <<SUMMARY
Input VCF: $INPUT
Source assembly: NCBI build 36 / hg18
Target assembly: GRCh38
Source PASS biallelic SNPs: $SOURCE_N
CrossMap unmapped records: $UNMAPPED
GRCh38 REF mismatches removed: $MISMATCH_REMOVED
Final GRCh38 PASS biallelic SNPs: $FINAL_N
Final C>T or G>A SNPs: $CTGA_N
Full catalogue: $FINAL
C-to-U-relevant catalogue: $CTGA
SUMMARY

cat "$OUTDIR/qc/summary.txt"
echo "[$(date '+%F %T')] COMPLETE"
echo "Recommended integration argument:"
echo "--variant-catalogue-vcf $CTGA"
