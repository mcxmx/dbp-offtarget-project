# v1.3.2 T227 Data Audit

Status: **NOT RECOVERED / NOT EVALUABLE as a row-level decomposition dataset (2026-09-17)**.

## Official-source search

The official SAMPDI-3D paper is Li et al., Bioinformatics 37:3760-3765, DOI `10.1093/bioinformatics/btab567`. The publisher/PMC supplementary PDF was recovered and inspected. It contains method and benchmark descriptions but no row-level T227 table, no per-mutation experimental values, and no per-mutation SAMPDI-3D or FoldX prediction table.

The paper reports an aggregate T227 benchmark of 227 DNA single-base-pair substitutions and an aggregate Pearson correlation of approximately 0.42 for SAMPDI-3D and 0.17 for FoldX. These aggregate values are not sufficient to compute global, position, or within-position identity decomposition and are not imported as pseudo-rows.

The article is internally inconsistent about cohort size: the methods/supplement describe 18 TFs, while the results text describes 17 TFs. This discrepancy is recorded rather than resolved by assumption.

The historical SAMPDI-3D web endpoint referenced by the paper (`compbio.clemson.edu/SAMPDI-3D`) is no longer available at the checked URL. A later PNBACE paper states that its Additional Table S3 is based on T227, but that derivative source is not the original SAMPDI-3D distribution and was not substituted for the frozen T227 audit in this stage.

## Endpoint status

| requested item | status | reason |
|---|---|---|
| 227 row-level experimental values | NOT_RECOVERED | absent from official paper supplement and unavailable from the historical server |
| TF/complex/PDB/mutation join | NOT_EVALUABLE | no row-level table |
| repeated mutant identities per position | NOT_EVALUABLE | cannot determine grouping without rows |
| published SAMPDI-3D per-row predictions | NOT_RECOVERED | only aggregate PCC reported |
| published FoldX per-row predictions | NOT_RECOVERED | only aggregate PCC reported |
| frozen global/position/identity decomposition | NOT_EVALUABLE | no auditable rows or deterministic mapping |

No T227 Parquet table or decomposition result is generated. The T227 claim boundary remains dataset recovery, not model performance.
