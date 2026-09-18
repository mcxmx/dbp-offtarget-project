# v1.3 data audit

Audit date: 2026-09-17. Counts below were computed from the checked-out files; no quantities were inferred.

## Scope and status

The competition source is `data/raw/v1_1_competition/source_data_extended_data_fig3.xls` (10 sheets). The seven proteins used throughout v0.x-v1.2 are development-exposed. DBP023, DBP056, and DBP062 occur in the source workbook and in the raw design-PDB directory but were absent from prior v1.1/v1.2 prediction and metric tables found by this audit; they are frozen as `LOCKED_HOLDOUT` before v1.3 method results. The primary method summary uses DBP001/003/005/006/009/035; DBP048 remains `SENSITIVITY_ONLY` because its mapping ambiguity was documented before v1.3.

The workbook contains a normalized PE/FITC value described as the mean of two replicates. It does not contain the two underlying replicate values. Therefore the v1.3 tidy table preserves one row per published mean and does not manufacture replicate rows. PBM has two replicate-level source files for the seven GSE237017 proteins.

At the original v1.3 freeze, `DBP35opt` was not present in the checkout. A v1.3.1 data-availability amendment recovered its 42-mutation, 14-position experimental landscape from official Extended Data Fig. 9e. The public workbook has a single summary column and no replicate-valued columns; DBP35opt model prediction remains unavailable.

## Per-protein inventory

| protein | assay | target sequence | target source | positions | variants | replicates | structure | analysis status | development/holdout | DeepPBS | missing |
|---|---|---|---|---:|---:|---|---|---|---|---|---|
| DBP001 | competition | `TAGCAGGATGTGT` | metadata/v0_3_1 experimental_assay_target | 13 | 39 | two-replicate mean only; raw replicates unavailable | available_file | PRIMARY | DEVELOPMENT_EXPOSED | available | individual competition replicates; DBP35opt pair |
| DBP003 | competition | `TAGCAGGATGTGT` | metadata/v0_3_1 experimental_assay_target | 13 | 39 | two-replicate mean only; raw replicates unavailable | available_file | PRIMARY | DEVELOPMENT_EXPOSED | available | individual competition replicates; DBP35opt pair |
| DBP005 | competition | `GCAGATCTGCACATC` | metadata/v0_3_1 experimental_assay_target | 14 | 42 | two-replicate mean only; raw replicates unavailable | available_file | PRIMARY | DEVELOPMENT_EXPOSED | available | individual competition replicates; DBP35opt pair |
| DBP006 | competition | `GCAGATCTGCACATC` | metadata/v0_3_1 experimental_assay_target | 14 | 42 | two-replicate mean only; raw replicates unavailable | available_file | PRIMARY | DEVELOPMENT_EXPOSED | available | individual competition replicates; DBP35opt pair |
| DBP009 | competition | `GCAGATCTGCACATC` | metadata/v0_3_1 experimental_assay_target | 14 | 42 | two-replicate mean only; raw replicates unavailable | available_file | PRIMARY | DEVELOPMENT_EXPOSED | available | individual competition replicates; DBP35opt pair |
| DBP023 | competition | `GCAGATCTGCACATC` | bundled PDB DNA chain B first strand | 14 | 42 | two-replicate mean only; raw replicates unavailable | available_file | LOCKED_HOLDOUT | LOCKED_HOLDOUT | not available | individual competition replicates; DBP35opt pair; PBM |
| DBP035 | competition | `GCAGATCTGCACATC` | metadata/v0_3_1 experimental_assay_target | 14 | 42 | two-replicate mean only; raw replicates unavailable | available_file | PRIMARY | DEVELOPMENT_EXPOSED | available | individual competition replicates; DBP35opt pair |
| DBP048 | competition | `CGACACCTGACGCG` | metadata/v0_3_1 experimental_assay_target | 14 | 42 | two-replicate mean only; raw replicates unavailable | available_file | SENSITIVITY_ONLY | DEVELOPMENT_EXPOSED | available | individual competition replicates; DBP35opt pair |
| DBP056 | competition | `CGCTATCCAGAGCG` | bundled PDB DNA chain B first strand | 14 | 42 | two-replicate mean only; raw replicates unavailable | available_file | LOCKED_HOLDOUT | LOCKED_HOLDOUT | not available | individual competition replicates; DBP35opt pair; PBM |
| DBP062 | competition | `CGCGATGCTTCTCG` | bundled PDB DNA chain B first strand | 14 | 42 | two-replicate mean only; raw replicates unavailable | available_file | LOCKED_HOLDOUT | LOCKED_HOLDOUT | not available | individual competition replicates; DBP35opt pair; PBM |

## Existing development artifacts

- uPBM: `data/processed/v0_3_1/designed_dbp_upbm_rc_class_v0_3_1.parquet`, 57,344 protein-RC-class rows, seven proteins, two PBM replicates per protein.
- Competition: `data/raw/v1_1_competition/source_data_extended_data_fig3.xls`, parsed into `data/processed/v1_3_competition_mutations.parquet`.
- DeepPBS: frozen NPZ bundle and v1.1 per-mutation table for seven proteins; DBP048 is mapping sensitivity-only in prior reports.
- Conditionality: frozen v0.6 M3 protein- and target-shuffle tables exist under `results/v0_6_representation/`; v1.3 standardizes them without rerunning predictions.
- Existing position/local metrics: `results/v1_2/04_position_decomposition/`, `results/v1_2/07_assay_triangle/`, and `results/v1_2/09_local_conditionality/`; v1.3 recomputes metrics from the tidy table and frozen predictions without changing their definitions after contract freeze.
- No DBP023/056/062 PBM or model prediction was found; these proteins are not silently promoted to an accuracy benchmark.

## Audit limitations

Structure-file presence is reported separately from validated structure/assay alignment. The three locked holdouts have raw PDB files but no checked-out v0.x-v1.2 structure metadata rows. Official Extended Data Fig. 3 was re-downloaded in v1.3.1 and contains the same 414 published means (maximum reconstruction error 0), but no underlying replicate values; competition replicate reproducibility therefore remains unavailable.
