# v0.5 Failure Analysis

## Scope and freeze

This is a secondary diagnostic over the frozen v0.5 primary evaluation. The
primary result files are protected by
`results/v0_5/PRIMARY_RESULTS_FROZEN_MANIFEST.txt` and were not overwritten.
Candidate-level scores were obtained by an exact replay of the frozen
fold/seed/config training protocol because Phase 3 did not persist predictions.
The replay was validated against the frozen seed-level Spearman values; the
maximum absolute difference was `9.714e-17`.

The hard-case reference is exactly the existing v0.3.1 set of **1,515**
sequence-vs-experiment disagreement candidates. No new candidate threshold was
chosen after inspecting M0-M3. Resolved means seed-mean prediction percentile
at least 0.90 within protein, matching the existing v0.4.2 rule.

## Hard-case resolution

| model | eligible | resolved | unresolved | resolution_rate |
| --- | --- | --- | --- | --- |
| M0 | 1515 | 383 | 1132 | 0.2528 |
| M1 | 1515 | 246 | 1269 | 0.1624 |
| M1c | 1515 | 282 | 1233 | 0.1861 |
| M2 | 1515 | 241 | 1274 | 0.1591 |
| M3 | 1515 | 236 | 1279 | 0.1558 |

Per-protein counts are in `results/v0_5/hard_case_resolution_by_protein.csv`.
The denominator is consistent across the five complete models: `1,515`
eligible candidates per model.

## Unique M3 wins and common failures

| subset | count | fraction_of_reference | n_proteins | proteins_represented | experimental_E_score_median | target_kmer_overlap_median |
| --- | --- | --- | --- | --- | --- | --- |
| m1_fail_m3_success | 164 | 0.1083 | 6 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6 | 0.2513 | 0.1429 |
| m1c_fail_m3_success | 154 | 0.1017 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2435 | 0.1429 |
| m2_fail_m3_success | 118 | 0.0779 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2286 | 0.1250 |
| joint_controls_fail_m3_success | 99 | 0.0653 | 6 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6 | 0.2257 | 0.1250 |
| all_current_models_fail | 871 | 0.5749 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2464 | 0.1667 |

`joint_controls_fail_m3_success` contains **99** candidates from
`DBP1|DBP3|DBP35|DBP48|DBP5|DBP6`. `all_current_models_fail.csv` contains **871**
candidates where M0, M1, M1c, M2, and M3 all remain unresolved under the
predeclared triage rule. These are ranking discrepancies, not biological
binding-failure claims.

The individual row-level files are:

- `m1_fail_m3_success.csv`
- `m1c_fail_m3_success.csv`
- `m2_fail_m3_success.csv`
- `joint_controls_fail_m3_success.csv`
- `all_current_models_fail.csv`

## Target-similarity artifact audit

| model | median_vs_target_kmer | median_vs_experimental |
| --- | --- | --- |
| M0 | 0.2280 | 0.1844 |
| M1 | 0.1021 | 0.0928 |
| M1c | 0.1074 | 0.0699 |
| M2 | 0.1718 | 0.0805 |
| M3 | 0.0717 | 0.1234 |

The per-protein and bin-level tables are
`target_similarity_model_correlations.csv` and
`target_similarity_bin_performance.csv`. Target-similarity tertiles were fixed
per protein before reading outcome values. M3's median prediction/target-k-mer
correlation was `0.0717`, compared with
M2 `0.1718` and M1 `0.1021`.

The residualized M3 diagnostic is in
`results/v0_5/partial_spearman_target_diagnostic.csv` and has median residual
association `0.1003`. This is a secondary rank
residual diagnostic, not a replacement for the primary Spearman metric.

## Protein and target use

| diagnostic | median_prediction_correlation | median_experimental_spearman | retrained |
| --- | --- | --- | --- |
| shuffled protein | 1.0000 | 0.0629 | 0 |
| shuffled target | 0.9999 | 0.0534 | 0 |

The shuffled-protein and shuffled-target diagnostics use already trained M3
models and set `retrained=False`; no shuffle result was used for training,
checkpoint selection, or tuning. Median prediction correlation after replacing
protein P was `1.0000`. After replacing target T it was
`0.9999`. A near-one correlation means the frozen
model's output is largely insensitive to that input; a large change without
better experimental ranking means the input affects predictions but does not
generalize in the intended direction.

## Why M0 can be stronger

| model | test_median_spearman | train_median_spearman |
| --- | --- | --- |
| M0 | 0.1549 | 0.2459 |
| M1 | 0.0571 | 0.2937 |
| M1c | 0.0910 | 0.3236 |
| M2 | 0.0608 | 0.2504 |
| M3 | 0.0629 | 0.2372 |

M0 has no cross-protein protein representation to overfit. The protein-aware
heads are trained on only the designed train proteins in each LOCO fold, while
the held-out protein is a new biological unit. Their median training/test
diagnostic values are in `train_vs_test_performance.csv`; training values are
descriptive only and were not used for selection.

Pair sampling used `512` pairs
per training protein. The median fraction of the 8,192 candidates appearing in
at least one sampled pair was `11.7%`, and
the median number of covered experimental rank deciles was
`10.0/10`. This is a plausible data-efficiency
bottleneck, but this phase does not increase pair counts.

M0's shared-signal diagnostics are in
`results/v0_5/m0_shared_signal_correlations.csv`, including correlations with
GC fraction and the prior k-mer3 proxy.

## SimplePC versus current M1

The prior v0.4.1 SimplePC result and current v0.5 M1 are not matched
experiments. SimplePC was trained on natural PBM data and evaluated on designed
proteins externally, while current M1 is trained only on the designed
protein-cluster LOCO training proteins. They also differ in protein
representation, DNA features, objective, pair sampling, normalization, and
evaluation regime. Therefore the apparent `0.3616` versus current
`0.0558`
gap is a protocol comparison, not evidence that one implementation is
intrinsically better.

No bridge experiment was run: it would require a new training comparison and
could be mistaken for an optimization after the frozen primary result.

## Decision

**MODIFY**

M3 does not beat the capacity-matched M1c or the target-blind M2 in the primary paired result, and the hard-case analysis does not provide sufficient evidence for a stable target-conditioned advantage. This rejects the current minimal implementation as the next training target, not the entire P,T,D concept.

This gate is about the current frozen minimal model family and evidence needed
for the next experiment. It does not establish that all target-conditioned
architectures are impossible.
