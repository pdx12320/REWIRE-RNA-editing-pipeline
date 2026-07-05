# Failure log

This file lists observed failures that changed the final workflow. It is not a generic troubleshooting list; each item corresponds to a problem encountered during development.

| Stage | Symptom | Root cause | Fix retained in final workflow |
|---|---|---|---|
| GATK preprocessing | `SAMRecord.getReadGroup() is null` | BAM header lacked a valid read group | verify or add `@RG` before `MarkDuplicates` |
| REDItools2 environment | Python import errors | Python 2 environment lacked compatible `mpi4py`, `pysam`, `sortedcontainers`, `psutil` or `netifaces` | use a dedicated Python 2 environment and pass its interpreter explicitly |
| MPI execution | insufficient slot error | requested workers exceeded allocated slots | inspect allocation and leave controller/system headroom |
| Supplementary contigs | failure around `chrGL000009` | temporary filename parsing removed `.1`/`.2` version suffixes | strip only the compression suffix and preserve full contig identifiers |
| Public WGS | 19.2–26.3% mapping and only 118 consensus variants | selected runs were unsuitable and represented different BioSamples | reject this route instead of treating it as genome-wide validation |
| Catalogue input | `Exec format error` | file content was gzip-compressed despite a `.vcf` filename | detect format from file bytes rather than extension |
| Chain retrieval | download timeout | remote chain download was unreliable on the server | accept and validate a local `hg18ToHg38.over.chain.gz` file |
| CrossMap output | tabix reported unsorted positions | liftover output was not coordinate-sorted | run `bcftools sort` before indexing |
| Catalogue liftover | 22,761 REF mismatches | lifted coordinates did not guarantee target-reference allele compatibility | validate against the project GRCh38 FASTA and remove mismatches |
| Final integration | missing `candidate_depth/` | legacy run did not preserve the newer depth-output directory | use the completed treatment-specific table only for supported exact-catalogue filtering |
| Final integration | missing `all_replicates_depth_pass` | site matrix came from an older helper implementation | do not fill or rename a surrogate column; use a dedicated compatibility script |

## Detailed notes

### 1. GATK read-group failure

A missing read group prevents GATK from associating reads with a sample. The fix is applied before duplicate marking, not after the BAM has already been processed.

Validation:

```bash
samtools view -H sample.bam | grep '^@RG'
```

### 2. REDItools2 Python 2 compatibility

REDItools2 required an older Python ecosystem. The final pipeline therefore passes the exact environment interpreter rather than relying on whichever `python` appears first in `PATH`.

```bash
"$CONDA_PREFIX/bin/python"
```

The same MPI implementation must be used by `mpirun`, `mpicc` and the compiled `mpi4py` module.

### 3. Supplementary-contig filename parsing

GRCh38 supplementary contigs can contain version suffixes. Removing text after the final period changed the contig identity. The retained parser removes only `.gz` and splits the remaining temporary filename from the right.

```python
pieces = os.path.basename(little_file)[:-3].rsplit("#", 2)
```

### 4. Public-WGS failure

The central problem was not merely low variant count. The resource had insufficient evidence to interpret absence as a reference genotype. The branch was therefore removed rather than weakened into a permissive blacklist.

### 5. Catalogue compression detection

The source file demonstrated that extensions are not reliable format metadata. The conversion script inspects the initial bytes:

```text
gzip: 1f 8b 08
text VCF: ##fileformat
```

### 6. Liftover and REF validation

CrossMap answers where a coordinate maps, but a biologically valid VCF record also requires the REF allele to match the target assembly. The retained workflow performs both coordinate liftover and allele validation.

### 7. Legacy integration boundary

The final catalogue became available after the earlier RNA filtering run had completed. Because the newer depth field was absent, the workflow was split into:

- a strict route for complete current inputs;
- a compatibility route for exact catalogue filtering of the already completed treatment-specific table.

This prevented a software-compatibility problem from being hidden as biological evidence.

## Rule derived from the failure log

Every automated safeguard in the final catalogue workflow corresponds to an observed failure:

```text
content detection
reference consistency
contig-name preservation
coordinate sorting
REF validation
exact-allele matching
explicit compatibility mode
```
