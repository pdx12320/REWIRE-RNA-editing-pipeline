import ast
import gzip
import re
from collections import defaultdict

BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def open_text(path, mode="rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else open(path, mode, newline="")


def norm_chrom(chrom):
    chrom = chrom.strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def parse_substitutions(raw, ref):
    if not raw or raw == "-":
        return []
    out = []
    for token in re.split(r"[\s,;]+", raw.strip()):
        token = token.upper().replace(">", "")
        if len(token) == 2 and token[0] == ref and token[1] in BASE_INDEX and token[1] != ref:
            out.append((token[0], token[1]))
    return out


def parse_counts(raw):
    try:
        values = ast.literal_eval(raw)
        if isinstance(values, (list, tuple)) and len(values) >= 4:
            return [int(float(x)) for x in values[:4]]
    except Exception:
        pass
    nums = re.findall(r"-?\d+(?:\.\d+)?", raw or "")
    if len(nums) >= 4:
        return [int(float(x)) for x in nums[:4]]
    raise ValueError("Cannot parse base counts: {!r}".format(raw))


def parse_vep(path):
    strand_sets = defaultdict(set)
    with open_text(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if not line.strip():
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            vid = row.get("Uploaded_variation", "")
            allele = row.get("Allele", "")
            strand = row.get("STRAND", "")
            if vid and strand in {"1", "-1"}:
                strand_sets[(vid, allele)].add(int(strand))
                strand_sets[(vid, "")].add(int(strand))
    return {k: next(iter(v)) if len(v) == 1 else 0 for k, v in strand_sets.items()}


def parse_wgs_vcf(path):
    variants = set()
    if not path:
        return variants
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or (len(p) > 6 and p[6] not in {".", "PASS"}):
                continue
            for alt in p[4].split(","):
                variants.add((norm_chrom(p[0]), int(p[1]), p[3].upper(), alt.upper()))
    return variants


def load_depth(path):
    depth = {}
    with open_text(path) as fh:
        for line in fh:
            if line.strip():
                chrom, pos, value = line.rstrip("\n").split("\t")[:3]
                depth[(norm_chrom(chrom), int(pos))] = int(value)
    return depth


def variant_id(chrom, pos, ref, alt):
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", chrom)
    return "REDI_{}_{}_{}_{}".format(safe, pos, ref, alt)
