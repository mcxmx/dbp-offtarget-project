# v1.3.3 Natural-TF Structural Model Evaluation

## Status

Primary natural-TF DeepPBS challenge: **NOT EVALUABLE**.

The contract and structure/mapping selection were frozen in `2e29f3c`; the
unscored prediction file was frozen in `e44e29e`. Only after that boundary were
the canonical SaMBA labels joined. The official DeepPBS preprocessing did not
produce an NPZ for the fixed 2STT ETS1 NMR input, so no model-vs-label
correlation or pairwise score is reported.

## Structure coverage

| item | count |
|---|---:|
| SaMBA TFs audited | 7 |
| SaMBA sites audited | 12 |
| structure-eligible sites | 4 |
| primary source-integrity-valid sites | 2 |
| eligible TFs | 1 (ETS1) |
| canonical mutation keys in prediction freeze | 132 |
| primary canonical mutation keys after v1.3.2 source-integrity gate | 66 |
| generated DeepPBS mutation scores | 0 |

The four structure-eligible sites are ETS1_site_1 through ETS1_site_4 mapped to
PDB 2STT by the label-blind deterministic sequence rule. The remaining eight
sites were excluded before prediction because the structure/DNA mapping was
ambiguous or did not cover every mutated coordinate. The frozen v1.3.2 source
integrity gate then restricts primary labels to ETS1_site_1 and ETS1_site_3;
ETS1_site_2 and ETS1_site_4 remain raw/audited but are not primary endpoints.
This is a structural/source-coverage limit, not a negative model result.

## Frozen decomposition table

`results/v1_3/natural_tf_decomposition_per_tf.tsv` contains one row per
eligible TF/site. All primary metrics are `NA` with status
`NOT_EVALUABLE_PREDICTION_NOT_REPRODUCED`. The summary table
`results/v1_3/natural_tf_decomposition.tsv` retains the pre-registered chance
baselines (0 for rank correlations and 0.5 for pairwise accuracy) and marks
all summaries `STRUCTURAL_COVERAGE_LIMITED_NOT_EVALUABLE`.

## Context for the four-level comparison

The source data for Figure 4 keep the quantities separate:

- designed development DeepPBS: model performance, development-exposed;
- designed DBP023/056/062: prospective holdout model performance;
- SaMBA canonical spot split: technical within-assay repeatability, not
  biological replication;
- natural-TF DeepPBS: not evaluable after the frozen primary preprocessing.

SaMBA technical repeatability remains high (TF-level median global 0.947,
position 0.881, identity residual 0.823, pairwise 0.814), but it is a
different assay estimand from competition and cannot be called a competition
noise ceiling.

## Decision

The natural-TF challenge cannot distinguish designed-protein domain shift from
a general structure-model identity limitation. The v1.4 identity-model gate
therefore remains **CLOSED**. No alternate PDB/register/orientation was chosen,
no DeepPBS version was substituted, and no GNN was trained.
