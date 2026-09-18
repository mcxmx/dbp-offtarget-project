# v1.3.1 Public Data Recovery and Prospective Holdout

## Decision summary

Official publisher source data recovered all ten published competition means and the DBP35opt landscape. They did not recover individual competition replicate values: the relevant downloadable sheets contain one summary column only. Replicate-level noise ceilings therefore remain not evaluable and no values were reconstructed from means or images.

The DeepPBS predictions for DBP023, DBP056, and DBP062 were generated without competition labels and frozen in commit `ff68394`. Evaluation occurred only after that commit, under the unchanged v1.3 contract.

| endpoint | holdout median | interpretation |
|---|---:|---|
| global Spearman | 0.085 | Weak aggregate mutation-effect ranking. |
| position Spearman | 0.560 | Consistent positive position-sensitivity signal. |
| identity residual Spearman | 0.001 | No detectable within-position identity ranking. |
| identity pairwise accuracy | 0.548 | Near the fixed 0.5 chance baseline. |

This is directionally consistent with Scenario D: `position >> identity` prospectively reproduces across the three held-out proteins. The protein-level sample size is three, so this is not treated as a broad significance claim.

## DBP35opt

The recovered DBP35opt experimental landscape is globally similar to DBP35 (Spearman 0.736) and retains position-sensitivity structure (Spearman 0.776), while within-position identity conservation is weaker (residual Spearman 0.184). The standard DBP35 and DBP35opt assays used different target and competitor concentrations, so the observed shift combines protein perturbation with assay-context difference. The model-shift endpoint remains unavailable because no DBP35opt prediction exists.

## Model-development gate

The evidence does not yet separate model-specific identity failure from weak identity-level assay identifiability because competition replicate values remain unavailable. The v1.4 identity-aware GNN gate is therefore not met. The highest-priority next step is a replicate-resolved competition source or an independent replicate-resolved mutation benchmark.
