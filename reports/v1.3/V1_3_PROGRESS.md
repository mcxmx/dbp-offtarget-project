# v1.3 progress

## Completed

- Audited README, PROGRESS, v1.1/v1.2 reports, scripts, results, metadata, raw competition workbook, PBM, structure files, and cached DeepPBS/NA-MPNN artifacts.
- Frozen `V1_3_DATA_AUDIT.md` and root `V1_3_EVALUATION_CONTRACT.md` before v1.3 benchmark results.
- Parsed all 10 competition sheets into `data/processed/v1_3_competition_mutations.parquet` and a consensus copy. Published values are two-replicate means; individual replicate rows are unavailable.
- Recomputed global, position-sensitivity, and within-position residual identity metrics for frozen DeepPBS/PBM/NA-MPNN inputs. Generated deterministic permutation null distributions and empirical p-values, protein bootstrap summaries, leave-one-protein-out sensitivity, conditionality diagnostics, and Figure 2 source data/PNG.
- Standardized the existing v0.6 M3 protein- and target-shuffle caches without rerunning inference; these remain development-exposed historical diagnostics.

## Key audit numbers

- Competition: 10 proteins, 414 mutation rows, 138 position records.
- Development-exposed: 7; locked holdout: 3 (DBP023/056/062).
- Primary reported method set: DBP001/003/005/006/009/035. DBP048 remains sensitivity-only because its local mapping ambiguity predates v1.3.
- DBP35opt: not found in checked-out data/metadata/results; real perturbation endpoint remains NOT EVALUABLE.
- DeepPBS: global median 0.354; position median 0.406; within-identity residual median 0.114; pairwise median 0.506 (chance 0.5).
- NA-MPNN: global median 0.149; position median 0.321; within-identity residual median 0.159; pairwise median 0.548 (chance 0.5).
- PBM_experimental: global median 0.364; position median 0.196; within-identity residual median 0.094; pairwise median 0.524 (chance 0.5).

## Stable conclusions

- The benchmark now distinguishes global ranking, position sensitivity, and within-position identity without pooling mutation rows as independent proteins.
- In the primary six, DeepPBS position-sensitivity agreement (median rho 0.406) is stronger than residual identity agreement (median rho 0.114); pairwise identity accuracy is 0.506, effectively the 0.5 chance baseline at the protein-summary level.
- The seven historical proteins remain development-exposed; DBP023/056/062 are not yet method-evaluable despite having competition rows.
- Cross-assay and competition-replicate variance decomposition is limited by absent individual competition replicate values; PBM replicate information is available.

## Not yet concluded

- DeepPBS DBP035 is the sole primary protein with nominal within-position evidence in both residual rho (0.444, permutation p=0.027) and pairwise accuracy (0.667, p=0.029). It is development-exposed and these p-values are not multiplicity-adjusted, so this is a case-level lead rather than validation.
- No claim is made about DBP35 -> DBP35opt because the paired protein is absent.
- No new GNN/model development is justified before reviewing the frozen decomposition and holdout availability.

## Next single priority

Obtain an auditable frozen prediction for DBP023/056/062 without using their competition values for model or metric tuning; otherwise the benchmark cannot make a prospective holdout claim.
