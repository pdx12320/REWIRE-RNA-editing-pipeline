# ORCA RNA-editing evidence pipeline

ORCA is the dry-lab system used to screen transcriptome-wide C-to-U RNA-editing candidates from three treated and three control RNA-seq libraries.

```text
RNA-seq alignment and preprocessing
→ REDItools2 substitution calling
→ transcript-strand interpretation
→ treated/control comparison
→ exact-allele filtering against a GRCh38-harmonized 293T catalogue
```

![ORCA system overview](wiki/assets/figure1_orca_system_overview.svg)

## Final result

| Evidence layer | Sites |
|---|---:|
| Strand-consistent site matrix | 9,930 |
| Called in all three treated replicates | 4,778 |
| Treatment-specific before catalogue comparison | 3,349 |
| Exact 293T catalogue overlaps | 16 |
| Final catalogue-filtered screening candidates | 3,333 |

![ORCA evidence funnel](wiki/assets/figure3_evidence_funnel.svg)

The retained sites are screening candidates rather than confirmed off-targets. The frozen legacy table does not contain independent depth and base counts for control non-calls, and the 293T catalogue is external to the exact experimental batch.

## Repository map

| Resource | Purpose |
|---|---|
| [Dry-lab wiki](wiki/README.md) | Finished ORCA workflow, results, contribution and limitations |
| [DBTL record](dbtl/README.md) | Ten development cycles plus failure and decision logs |
| [Pipeline guide](pipeline/README.md) | Reproducible commands |
| [Catalogue provenance](pipeline/CATALOGUE_PROVENANCE.md) | Source, liftover and quality control |
| [Outputs](pipeline/OUTPUTS.md) | Expected files |
| [Troubleshooting](pipeline/TROUBLESHOOTING.md) | Observed errors and fixes |
| [Result summary](results/README.md) | Frozen counts |

## DBTL cycles

```text
1  RNA evidence
2  public WGS test
3  catalogue harmonization
4  final integration
5  GATK read groups
6  REDItools2 environment
7  contigs and coverage
8  control-depth evidence
9  legacy compatibility
10 audit and reporting
```

Large sequencing and intermediate files remain outside GitHub. The repository retains scripts, environment definitions, quality-control rules, frozen counts and decision records.

> Repository note: the GitHub URL retains the earlier repository name for continuity, but the system described in the documentation is ORCA.
