# REDItools2 resume and recovery

The repository runner is designed to recover from two upstream REDItools2 failure modes observed with GRCh38:

1. supplementary contig identifiers such as `GL000009.2` and `KI270750.1` are shortened incorrectly during the final temporary-file sort;
2. `files.txt` is missing even though every interval file has already been produced.

These errors occur after the expensive MPI analysis can already be complete. Rerunning the sample from the beginning is therefore unnecessary when all expected interval files are present.

## What the runner now does

`scripts/run_reditools_all_samples.sh` now:

- passes a temporary directory with a trailing slash;
- checks for an already valid compressed table and tabix index;
- validates legacy coverage output before reusing it;
- attempts to rebuild `files.txt` from completed interval files before starting MPI again;
- preserves exact contig identifiers, including `.1`, `.2`, or any other version suffix;
- compares the number of recovered interval files with the expected count stored in `intervals.txt`;
- merges completed temporary output even when the upstream process returned a non-zero status during its final sort;
- writes the final result atomically and verifies its gzip stream, header, and tabix index;
- archives incomplete temporary directories instead of silently mixing old and new interval files.

## Recovery helper

The file list is rebuilt with:

```bash
python3 scripts/rebuild_reditools_file_list.py \
  --temp-dir "$PROJECT/reditools/tmp/CU517_GC_T2" \
  --fai "${REF}.fai"
```

The helper removes only the terminal `.gz` extension from a chunk filename. For example:

```text
GL000009.2#100#50000.gz
```

is parsed as:

```text
contig = GL000009.2
start  = 100
end    = 50000
```

The `.2` suffix is retained because it is part of the GRCh38 contig identifier.

## Resume command

After pulling the latest repository version:

```bash
cd /data/ydx/igem/REWIRE-RNA-editing-pipeline
git pull
conda activate reditools2_py2

PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
REDITOOLS=/data/ydx/igem/REDItools2

nohup bash scripts/run_reditools_all_samples.sh \
  "$PROJECT" \
  "$REF" \
  "$REDITOOLS" \
  "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > "$PROJECT/logs/reditools_resume.log" 2>&1 &
```

If a sample has complete temporary interval output but no final table, the script rebuilds `files.txt` and merges it. It does not rerun the expensive MPI stage.

## Validation

A completed sample must pass all of the following:

```bash
gzip -t "$PROJECT/reditools/tables/CU517_GC_T2.txt.gz"
tabix -l "$PROJECT/reditools/tables/CU517_GC_T2.txt.gz" | head
zcat "$PROJECT/reditools/tables/CU517_GC_T2.txt.gz" | head -n 3
```

The first row must contain the standard 14-column REDItools2 header.
