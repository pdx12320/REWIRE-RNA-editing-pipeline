# Decision log

This file records the major methodological choices, alternatives considered and evidence used to select the final implementation.

## D1 — Analyse replicates independently

**Decision:** process each of the three treated and three control libraries separately.

**Alternative considered:** pool treated reads and pool control reads before calling.

**Reason:** separate calls retain replicate reproducibility and expose sample-specific failures. Pooling could allow one high-depth replicate to dominate the result.

**Consequence:** final tables preserve replicate-level coverage, alternate-read count and editing rate.

---

## D2 — Use transcript orientation rather than genomic substitution alone

**Decision:** interpret genomic C>T on positive-strand transcripts and genomic G>A on negative-strand transcripts as transcript-level C-to-U candidates.

**Alternative considered:** retain only genomic C>T calls.

**Reason:** a C-to-U change on a negative-strand transcript appears as G>A in genomic coordinates. Ignoring strand would systematically discard valid negative-strand candidates.

**Consequence:** VEP strand annotation became a required integration input.

---

## D3 — Keep a stringent REDItools2 edited-read threshold

**Decision:** retain the original `-me 20` edited-read threshold.

**Alternative considered:** lower the threshold to increase sensitivity.

**Reason:** the first objective was to generate a strongly supported screening set. Lower thresholds would increase sensitivity but also increase the burden of distinguishing sequencing noise and low-level background.

**Boundary:** the resulting call set may miss low-frequency editing.

---

## D4 — Do not interpret a missing control call as zero

**Decision:** describe control non-calls as missing call evidence, not confirmed zero editing.

**Alternative considered:** assign edit rate 0 to every control non-call.

**Reason:** REDItools2 can omit a site that has lower-level alternate reads but does not reach the edited-read threshold.

**Consequence:** the strict workflow includes direct candidate-site depth measurement. The frozen legacy result is labelled as a screening set because the independent control-depth outputs were unavailable at final integration.

---

## D5 — Reject the public-WGS blacklist

**Decision:** remove the three-public-run WGS route from the final workflow.

**Alternative considered:** retain the 118-site two-of-three consensus as a genomic blacklist.

**Evidence:**

```text
mapping rate: 19.2–26.3%
strict per-sample calls: 137–230
2-of-3 consensus: 118 variants
```

**Reason:** absence from such a sparse call set would reflect missing genomic evidence, not a reference genotype.

**Consequence:** the project moved to the HEK293 Genome Project `293T_CG` catalogue.

---

## D6 — Use an external 293T catalogue with explicit scope limits

**Decision:** use the database-released `293T_CG` call set as external exclusion evidence.

**Alternative considered:** claim it as matched WGS.

**Reason:** the catalogue is genome-scale and biologically relevant to 293T, but it was not generated from the exact CU5.17 experimental batch.

**Consequence:** exact overlaps are treated as plausible genomic variants; absence is not treated as proof that the experimental cells are variant-free.

---

## D7 — Convert hg18 to GRCh38 before any comparison

**Decision:** perform formal coordinate liftover and REF validation.

**Alternative considered:** compare positions directly or use a coordinate-only mapping table.

**Reason:** the RNA branch uses GRCh38 while the source catalogue uses build36/hg18. Direct comparison would be invalid.

**Consequence:** CrossMap, the UCSC chain and project FASTA validation are mandatory.

---

## D8 — Validate REF after liftover

**Decision:** remove records whose lifted REF allele does not match the project GRCh38 FASTA.

**Alternative considered:** retain every successfully mapped coordinate.

**Evidence:** 22,761 mapped records had target-reference mismatches.

**Reason:** coordinate mapping alone does not guarantee a valid target-assembly VCF allele.

**Consequence:** `bcftools norm -f ... -c x` is part of the frozen catalogue contract.

---

## D9 — Match exact alleles, not positions

**Decision:** use `CHROM:POS:REF:ALT` identity.

**Alternative considered:** mark every RNA candidate sharing a catalogue coordinate.

**Reason:** the same coordinate can contain different alternate alleles. Coordinate-only matching would over-filter.

**Consequence:** chromosome naming is normalized, while REF and ALT remain required components of the key.

---

## D10 — Preserve excluded records

**Decision:** write catalogue-overlapping candidates to a separate output table.

**Alternative considered:** delete them during filtering.

**Reason:** an auditable pipeline should show what was removed and why.

**Consequence:** the final output includes both retained and excluded tables plus a count summary.

---

## D11 — Add a compatibility route instead of rewriting history

**Decision:** add `filter_existing_treatment_specific_by_293T.py` for the completed legacy table.

**Alternative considered:** fill the missing `all_replicates_depth_pass` field, rename another column or rerun the entire pipeline without verified original inputs.

**Reason:** the available legacy table supported exact catalogue comparison but did not support reconstruction of independent control-depth evidence.

**Consequence:** the final 3,333 records are explicitly called screening candidates.

---

## D12 — Separate the public-facing wiki from the development record

**Decision:** keep the wiki focused on the finished dry-lab workflow, results, contribution and limitations; store DBTL history in this folder.

**Alternative considered:** place every failure and implementation detail on the wiki page.

**Reason:** readers need the biological question and result first, while reproducibility requires a deeper record that remains accessible without interrupting the main narrative.

**Consequence:**

```text
wiki/README.md        finished project story
dbtl/                 iterative development record
pipeline/             executable implementation
```
