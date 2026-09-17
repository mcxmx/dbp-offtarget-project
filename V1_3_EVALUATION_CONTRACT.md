# V1.3 Evaluation Contract

Status: **FROZEN 2026-09-17**. This contract is written before v1.3 results. Metric definitions may change only for a demonstrated implementation bug; a result that is inconvenient is not a bug.

## Analysis units and independence

The primary statistical unit is the protein. Mutation rows are repeated observations within a protein landscape, not independent biological replicates. All primary summaries report per-protein values, the number of proteins, median, mean as an auxiliary summary, protein-level bootstrap 95% intervals, and leave-one-protein-out values where implemented. Development-exposed proteins are not independent validation.

The seven GSE237017 proteins and all v1.1/v1.2 model outputs are `development_exposed`. DBP023, DBP056, and DBP062 are frozen `LOCKED_HOLDOUT` because the source workbook contains them but the audit found no prior prediction/metric artifact using them. They cannot be used for tuning. This status is frozen before reading v1.3 benchmark results.

The primary method summary contains DBP001, DBP003, DBP005, DBP006, DBP009, and DBP035. DBP048 remains `SENSITIVITY_ONLY` because v1.1/v1.2 had already documented an ambiguous local mapping; this exclusion was inherited before v1.3 outcomes and is not result-driven. An all-seven descriptive summary is retained and explicitly labeled.

## Measurement direction and ties

Competition raw measurement is normalized PE/FITC. The frozen monotone effect is `normalized_effect = -raw_measurement`; larger values mean stronger competition. The raw value remains in the tidy table. PBM and competition are different assay estimands and are never treated as interchangeable ground truth. Spearman uses average ranks (the scipy default). Pairwise identity accuracy gives a comparable pair score of 1 for a correct strict order, 0 for an incorrect strict order, and 0.5 for a prediction tie; experimental ties are excluded from the denominator. Its chance baseline is 0.5.

## A. Global mutation-effect endpoint

For every protein and frozen method, calculate Spearman and Pearson over all available mutation rows for that protein. Report the per-protein values. A descriptive rank-based 95% interval is obtained by deterministic position-stratified bootstrap of mutation rows within that protein; it is not a biological-replicate interval. The primary across-protein interval is a protein bootstrap (2,000 resamples, seed 1301) of the per-protein statistic. No pooled mutation-row p-value is primary.

## B. Position sensitivity endpoint

For protein p and position j, the primary position sensitivity is `mean_abs_effect(j) = mean_b(abs(effect(j,b)))` over the three substitutions. This definition is fixed because competition direction is an assay-specific signed signal and absolute magnitude directly represents sensitivity without choosing a deleterious sign threshold. A WT-relative "mean deleterious effect" is not identifiable from this workbook because it provides normalized competitor signal rather than a WT-relative signed change. `mean_signed_effect` and `max_abs_effect` are retained as secondary diagnostics. Predicted values use the same transformation. The endpoint is Spearman between experimental and predicted position summaries, per protein, with protein-level bootstrap intervals.

## C. Within-position nucleotide identity endpoint

For each position, residualize both vectors by their within-position mean: `delta(j,b) = effect(j,b) - mean_b(effect(j,b))`. Primary metrics are pooled within-position residual Spearman and Pearson within each protein, plus mean pairwise accuracy across positions. Identity pair comparisons are only within a position. Position-centering removes position sensitivity; no global centering is used. Chance is 0 for residual correlation under a symmetric random-order null and 0.5 for pairwise accuracy. The reported Kendall-style concordance is the equivalent transform `2 * pairwise_accuracy - 1`, with chance 0; it is not treated as an additional independent endpoint. Within-position base-label permutation shuffles predicted mutant labels within each position (2,000 permutations, seed 1301); this is the primary identity null.

## D. Protein conditionality

For cached global landscapes, report pairwise protein Spearman, mean absolute prediction shift between proteins, and per-7-mer variance across protein identities. The landscape conditionality effect size is the mean across-7-mer between-protein variance divided by total variance across all protein/7-mer scores; it is descriptive and is not an accuracy claim. For the existing frozen v0.6 M3 protein- and target-shuffle predictions, preserve prediction correlation, mean absolute score shift, and the historical normalized rank-shift definition `1 - prediction_correlation`. These are development-exposed diagnostics, not a new holdout evaluation. Protein-label permutation is at the protein unit. If a future method has no target-conditioned artifact, target shuffle is reported unavailable rather than simulated.

The real perturbation endpoint requires matched DBP35 and DBP35opt experimental and predicted landscapes. The audit found no DBP35opt source, so this endpoint is pre-specified but currently `NOT EVALUABLE`, not replaced by a synthetic proxy.

## Nulls and uncertainty

Position-label permutation shuffles predicted position summaries across positions within protein. Within-position base-label permutation shuffles predicted mutant identity within each position. Protein-label permutation is applied only at the protein unit. Null distributions and observed values are written to source tables. Correlation endpoints use two-sided empirical permutation p-values around their chance value of zero; pairwise accuracy uses the upper tail around chance 0.5. The add-one correction is used. Leave-one-protein-out sensitivity is descriptive; no protein is selected after seeing its result.

## Method inclusion

The first benchmark includes frozen DeepPBS, frozen NA-MPNN landscape-derived scores when the local register has coverage, experimental PBM as an assay comparator, and position-only summaries. No new neural model is trained in v1.3. Methods lacking a complete, auditable prediction are marked unavailable in the method reproducibility record and are not backfilled.

## Figure 2 data contract

Figure 2 source data must contain one row per method/protein with global Spearman, position-sensitivity Spearman, within-position residual Spearman, within-position pairwise accuracy, mutation count, holdout/development status, and the fixed chance baselines. The figure is generated only from this source table; no hand-edited numbers are permitted.
