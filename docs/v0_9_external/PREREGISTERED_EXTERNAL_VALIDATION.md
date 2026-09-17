# Preregistered External Validation

**Status:** FROZEN BEFORE LABEL ACCESS. No independent quantitative labels have been read. The current project state is `WAIT FOR DATA; CONTACT AUTHORS`.

## Hypothesis

Under strict unseen-protein evaluation, frozen sequence-based specificity models will show lower transfer performance on an independent de novo DBP cohort than on matched experimental replicate references; structure-derived methods, when input-compatible, are evaluated descriptively rather than assumed to rescue the gap.

## Cohort inclusion

The cohort must contain de novo designed DNA-binding proteins, exact experimental construct sequences, intended targets, quantitative specificity measurements, multiple independent designs (preferably more than seven), and a Tier 1 or Tier 2 assay. Tier 3 assays are eligible only as orthogonal validation. The assay schema, design IDs, randomized region and protein-level coverage must be inspectable without opening labels.

## Exclusion

Exclude development-exposed GSE237017/Glasscock data, missing construct or target mapping, incompatible assay, ambiguous register, unresolvable strand convention, insufficient protein-level coverage, or files whose schema cannot be audited before label access. Do not silently drop an eligible protein after labels are visible.

## Frozen methods and scoring

Allowed methods and checkpoints are in `metadata/v0_9_external/external_validation_config.yaml`: sequence k-mer, corrected v0.6 M0 and M3, SimplePC only if assay/protocol compatible, DeepPBS and NA-MPNN only with compatible complex structures, and AF3/contact probability only if its predefined subset and compute contract are met. v0.6 seeds are aggregated by the fixed mean of seeds 17, 29 and 43. Structure PWM/PPM outputs use the fixed log-probability proxy inherited from v0.8.

Canonical reverse-complement classes are used. DNA register and strand convention are declared from assay metadata before labels; an ambiguous register excludes rather than triggers an alignment search. Tier-specific endpoints are exactly those in `ASSAY_SPECIFIC_VALIDATION_RULES.md`.

## Analysis

Primary summaries are median per-protein endpoint values and complete per-protein tables. Secondary analyses include predeclared top-k/NDCG or discrimination metrics, coverage, failures and exploratory conditionality diagnostics where meaningful. Protein is the independent statistical unit. Use paired exact/permutation summaries and effect sizes; never use candidate count as the population N.

## Missing data and failures

Report every eligible protein, missing fields, execution errors and coverage. Preserve failed outputs and logs. There is no automatic retry with another scoring rule, register, checkpoint or hyperparameter.

## No post-hoc selection

After label access, do not change model set, checkpoint, seed aggregation, register, scoring function, primary metric, exclusions or protein denominator. Any additional analysis is marked `POST-HOC / EXPLORATORY`.

## Confirmatory status

The single frozen pipeline is confirmatory. External validation has not started because no cohort has yet passed the metadata firewall. The accompanying manifest records `label_access: false`.

The runner is deliberately fail-closed while the cohort schema and method adapters are unverified: it refuses to write primary results rather than emitting a placeholder or silently changing an estimand.
