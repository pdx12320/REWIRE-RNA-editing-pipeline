# Postprocessing commands

After six REDItools2 tables are complete:

```bash
PROJECT=/data/ydx/igem/CU5.17_EGFP_GC_paper
MANIFEST=config/samples.tsv

mkdir -p "$PROJECT/vcf"
python3 scripts/reditools_union_to_vcf.py \
  --manifest "$MANIFEST" \
  --reditools-dir "$PROJECT/reditools/tables" \
  --vcf "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.vcf" \
  --bed "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed"

python3 scripts/run_vep_annotation.py \
  --project "$PROJECT" \
  --cache /path/to/vep_cache

bash scripts/build_candidate_depth_tables.sh \
  "$PROJECT" \
  "$PROJECT/vcf/CU5.17_EGFP_GC.REDItools_union.bed" \
  "$MANIFEST"
```

Then run `filter_c_to_u_and_compare.py` with the manifest, REDItools2 table directory, VEP table, candidate-depth directory, final output directory, and optional matched WGS VCF.
