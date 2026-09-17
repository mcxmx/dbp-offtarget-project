# v1.3 method reproducibility

This record freezes method inclusion before interpreting Figure 2. A method is included in the local decomposition only when a cached prediction can be joined to the competition key `(protein, position, wt_base, mutant_base)` without refitting or selecting a register from the v1.3 outcome.

| method | local input | status | reason |
|---|---|---|---|
| DeepPBS | `results/v1_1/05_assay_aligned_eval/deeppbs_per_mutation_predictions.tsv` | INCLUDED | Frozen native PWM delta-logP; all seven development-exposed proteins, fixed v1.1 register. |
| NA-MPNN | `results/v1_1/05_assay_aligned_eval/nampnn_local_per_mutation_predictions.tsv` | INCLUDED_WITH_COVERAGE_LIMIT | Frozen global-landscape-derived local scores; only the 7 positions covered by the fixed mapped window (21 mutations/protein). No new inference. |
| PBM experimental comparator | `results/v1_1/06_assay_agreement/pbm_competition_per_mutation.tsv` | INCLUDED_AS_ASSAY_COMPARATOR | Experimental PBM delta from the frozen 7-mer register; not a computational method and not interchangeable ground truth. |
| M0/M1/M1c/M2/M3 | historical v0.5-v0.6 global 7-mer artifacts | CONDITIONALITY_CACHE_ONLY | Available artifacts do not expose a frozen competition-key prediction table for Figure 2. Existing M3 protein/target shuffle tables are standardized in `results/v1_3/cached_shuffle_diagnostics.tsv` without new inference. |
| historical sequence baseline/SimplePC | v0.3-v0.4 benchmark artifacts | NOT_INCLUDED_IN_LOCAL_FIGURE_2 | Sequence/PBM benchmark outputs are not a frozen per-mutation prediction table for this competition assay and cannot be silently treated as a local model. |
| DBP023/DBP056/DBP062 DeepPBS | `results/v1_3/holdout_predictions_unscored.tsv` | PROSPECTIVE_HOLDOUT_EVALUATED_AFTER_FREEZE | Prospective inference used official commit `8bfb211`, the existing five-checkpoint ensemble, design PDB chains A/B, and a PDB-order mapping fixed before label access. Stage 1 was committed as `ff68394` before labels were joined; Stage 2 results are in `results/v1_3/holdout_decomposition.tsv`. |
| DBP35 -> DBP35opt experiment | official Extended Data Fig. 9e | EXPERIMENT_RECOVERED_CONTEXT_CONFOUNDED | The official numeric landscape is available, but DBP35opt used 20 nM target/160 nM competitor versus the standard 1 uM/8 uM condition. No DBP35opt model prediction is available, so predicted landscape shift remains not evaluable. |

The v1.3 benchmark itself used cached predictions without rerunning external methods. The v1.3.1 prospective holdout is the explicitly documented exception: it ran the same official DeepPBS checkout and frozen pipeline for three previously unscored structures, without holdout labels or parameter selection. Prediction generation and evaluation are separated by freeze commit `ff68394`; old v0.x-v1.2 outputs were not overwritten.
