# Current REDItools2 run

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
REF=/data/ydx/igem/GRCh38.primary_assembly.genome.fa
REDITOOLS=/data/ydx/igem/REDItools2

conda activate reditools2_py2

nohup bash scripts/run_reditools_all_samples.sh \
  "$PROJECT" "$REF" "$REDITOOLS" \
  "$CONDA_PREFIX/bin/python" \
  30 8 8 \
  > reditools.log 2>&1 &
```

The last three numbers are MPI processes, concurrent coverage jobs, and compression threads. Monitor with `tail -f reditools.log`.
