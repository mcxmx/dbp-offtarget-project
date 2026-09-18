# Supplementary Material Plan

## Tables

1. **Supplementary Table 1.** All designed proteins, cohort assignment, assay status, mapping status, and development/holdout status.
2. **Supplementary Table 2.** Per-protein DeepPBS global, position, identity-residual, pairwise, and permutation metrics.
3. **Supplementary Table 3.** Locked holdout predictions and decomposition for DBP023, DBP056, and DBP062, including the pre-reveal freeze record.
4. **Supplementary Table 4.** Canonical SaMBA technical spot-split repeatability by TF/site and the across-TF summary.
5. **Supplementary Table 5.** 2STT structure audit, rescue routes, process exits, expected NPZ path, and final reproducibility status.
6. **Supplementary Table 6.** Per-protein experimental and predicted variance decomposition and position-only counterfactual values.

## Figures

1. Per-protein permutation distributions for global, position, residual identity, and pairwise metrics.
2. Per-protein position-only counterfactual scatterplots and score deltas.
3. DBP35/DBP35opt context-confounded landscape comparison.
4. Natural-TF technical rescue audit and MODEL 1 provenance.
5. Leave-one-protein-out summaries for development-set sensitivity.

## Reporting rules

- SaMBA is labelled **technical repeatability**, never biological replicate or competition noise ceiling.
- The natural-TF result remains **NOT_EVALUABLE**; any later MODEL 1 output is **POST-REVEAL TECHNICAL SENSITIVITY**.
- Variance fractions describe landscape allocation, not the fraction of predictive performance caused by a feature.
- DBP35opt remains a context-confounded case study and is not interpreted as a matched causal mutation experiment.
- No new model training is part of the v1.3.4 supplement.
