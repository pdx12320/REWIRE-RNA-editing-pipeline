# Cycle 2 — Test public WGS as a genomic blacklist

## Design

### Question

Could public 293T whole-genome sequencing be used to remove RNA mismatches that were actually genomic variants?

RNA-seq alone cannot reliably distinguish an editing event from a DNA variant. The first solution considered was to reconstruct a genomic blacklist from public WGS data and remove exact RNA alleles found in that blacklist.

### Initial plan

Three public WGS runs were evaluated:

```text
SRR37832939
SRR37832940
SRR37832941
```

The runs represented different BioSamples rather than technical replicates from one experimental cell batch. Therefore, the planned integration rule was a two-of-three consensus rather than merging all reads into one sample.

The intended workflow was:

```text
public WGS FASTQ
→ GRCh38 alignment
→ BAM quality control
→ per-sample variant calling
→ strict PASS biallelic SNP filtering
→ two-of-three exact-allele consensus
→ genomic blacklist for RNA candidates
```

## Build

Each WGS run was handled independently to avoid treating genetically distinct BioSamples as replicates of one genome.

The evaluation focused on:

- fraction of reads mapping to GRCh38;
- usable depth after alignment;
- number of filtered variants per run;
- size of the exact-allele consensus;
- whether the resulting call set was plausible for a human genome.

A useful blacklist would need broad genome-wide coverage and millions of callable sites or variants. A small number of calls would indicate that absence from the WGS output was mostly missing evidence rather than evidence of a reference genotype.

## Test

### Mapping performance

Only 19.2–26.3% of reads mapped to GRCh38 across the three runs.

This was far below what would be expected for suitable human WGS and indicated that the selected public data were not a reliable matched genomic resource for this project.

### Variant yield

Strict per-sample filtering retained only:

```text
137–230 variants per sample
```

The two-of-three consensus contained:

```text
118 variants
```

For comparison with the later catalogue workflow, the HEK293 Genome Project `293T_CG` call set retained 2,885,725 GRCh38 PASS biallelic SNPs after conversion and reference validation. The public-WGS result was therefore several orders of magnitude too sparse.

### Interpretation

A candidate RNA allele absent from a 118-site blacklist cannot be considered genomically excluded. Most of the genome was simply uninformative.

The low yield also meant that retaining the route would create a misleading appearance of genomic validation without providing meaningful coverage.

## Learn

### Decision

```text
reject the public-SRA WGS blacklist
→ do not use the 118-site consensus in the final analysis
→ search for a database-released 293T genomic catalogue
```

### Why the route was removed

1. The mapping rate was inadequate.
2. The three runs were not matched to the exact CU5.17 experimental batch.
3. They represented different BioSamples.
4. Variant counts were implausibly small for genome-wide exclusion.
5. A negative result from this resource would mainly reflect missing coverage.

### Reporting rule

The failed WGS branch is documented because it changed the project design, but it is not presented as part of the final validation evidence.

The final wiki states the decision briefly; the full diagnostic history remains here so that later users do not repeat the same route without checking sample identity and genome-wide callability.

## Output of this cycle

No WGS blacklist from these three runs was retained.

The next cycle replaced this approach with the database-released `293T_CG` VCF from the HEK293 Genome Project.
