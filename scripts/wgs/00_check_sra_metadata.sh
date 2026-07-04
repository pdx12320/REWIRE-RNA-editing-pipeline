#!/usr/bin/env bash
set -Eeuo pipefail

RUNS_FILE="${1:-config/wgs_runs.tsv}"
OUT="${2:-sra_metadata.tsv}"

command -v curl >/dev/null 2>&1 || { echo "ERROR: curl not found" >&2; exit 1; }
[[ -s "$RUNS_FILE" ]] || { echo "ERROR: run manifest not found: $RUNS_FILE" >&2; exit 1; }

FIELDS='run_accession,study_accession,sample_accession,experiment_accession,scientific_name,library_strategy,library_source,library_layout,instrument_platform,instrument_model,base_count,read_count'
first=1
: > "$OUT"

while IFS=$'\t' read -r run; do
    [[ -z "$run" || "$run" == "run" || "$run" == \#* ]] && continue
    url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=${FIELDS}&format=tsv"
    tmp="$(mktemp)"
    if ! curl -fsSL --retry 4 --retry-delay 5 "$url" -o "$tmp"; then
        rm -f "$tmp"
        echo "ERROR: failed to retrieve ENA metadata for $run" >&2
        exit 1
    fi
    if [[ $(wc -l < "$tmp") -lt 2 ]]; then
        cat "$tmp" >&2
        rm -f "$tmp"
        echo "ERROR: no metadata returned for $run" >&2
        exit 1
    fi
    if [[ $first -eq 1 ]]; then
        cat "$tmp" >> "$OUT"
        first=0
    else
        tail -n +2 "$tmp" >> "$OUT"
    fi
    rm -f "$tmp"
done < "$RUNS_FILE"

echo "Metadata written to: $OUT"
column -t -s $'\t' "$OUT" 2>/dev/null || cat "$OUT"

echo
echo "Interpretation:"
python3 - "$OUT" <<'PY'
import csv, sys
path=sys.argv[1]
with open(path, newline='') as fh:
    rows=list(csv.DictReader(fh, delimiter='\t'))
if not rows:
    raise SystemExit('No metadata rows found')
strategies={r.get('library_strategy','').upper() for r in rows}
layouts={r.get('library_layout','').upper() for r in rows}
samples={r.get('sample_accession','') for r in rows if r.get('sample_accession','')}
organisms={r.get('scientific_name','') for r in rows}
print('Runs:', ', '.join(r['run_accession'] for r in rows))
print('Library strategy:', ', '.join(sorted(strategies)))
print('Layout:', ', '.join(sorted(layouts)))
print('BioSamples:', ', '.join(sorted(samples)) or 'not reported')
print('Organism:', ', '.join(sorted(organisms)))
if strategies != {'WGS'}:
    print('WARNING: not all runs are labelled WGS.')
if layouts != {'PAIRED'}:
    print('WARNING: not all runs are paired-end.')
if len(samples) == 1 and len(rows) == 3:
    print('Recommended mode: merge (three runs belong to one BioSample).')
else:
    print('Recommended mode: consensus (runs are different BioSamples or metadata are incomplete).')
PY
