# v0.4.2 Baseline Gap Analysis

Date: 2026-09-02

## Main ranking results

| Method | Designed macro-median Spearman | Natural-test macro-median Spearman | Coverage |
|---|---:|---:|---:|
| Best sequence-only k-mer3 proxy | 0.232 | NA | 7/7 |
| SimpleProteinConditionalBaseline | 0.362 | 0.301 | 7/7 |
| FrozenPLMProteinConditionalBaseline | 0.153 | 0.316 | 7/7 / natural test |
| DeepPBS | NA | NA | 0/7 |
| NA-MPNN diagnostic | partial | NA | 2/7 |
| Empirical replicate reference | 0.591 | NA | 7/7 reference |

The replicate value is an empirical assay reproducibility reference, not a
strict theoretical upper bound. Scores are compared per protein before
macro-summary; raw scores are not pooled across proteins.

## Failure diagnostics

- The pre-registered v0.3.1 disagreement set contains 1,515 RC-class candidates.
- Resolution counts under the v0.4.2 top-10% rule are recorded in
  `tables/disagreement_resolution_v0_4_2.csv`; missing predictions are
  `not_evaluable`, never zero-filled.
- The complete-core common-hard set contains 263 rows:
  experimental percentile >= 0.95 and sequence k-mer3, SimplePC, and FrozenPLM
  prediction percentiles <= 0.50.
- Common-hard rows are distributed across all seven DBPs:
  DBP1=40, DBP3=42, DBP35=30, DBP48=72, DBP5=14, DBP6=30, DBP9=35.
- The lowest FrozenPLM per-protein result is DBP48 at
  -0.437 Spearman. DBP48 is also negative
  (-0.437).
- Motif-distance stratification is descriptive and does not prove epistasis.
  The output is `tables/performance_by_motif_distance_v0_4_2.csv`.

## Interpretation

FrozenPLM provides a stronger protein representation than the composition
baseline on the natural test split, but its designed external median is below
both the SimplePC result and the best sequence-only k-mer3 result. This is a
real baseline result, not evidence for a proposed architecture. The current
gap analysis is insufficient for a final model-development GO because DeepPBS
was not executed in a Linux-compatible runtime and the assay-aligned natural
construct benchmark remains empty.
