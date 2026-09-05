# v0.5 Primary Evaluation Results

## Experimental status

This report is the frozen v0.5 primary evidence artifact. It reports the complete four-fold protein-cluster LOCO evaluation and the three-fold combined-component sensitivity evaluation using the fixed seeds 17, 29, and 43.

`protein_cluster_loco_fold_1` was previously used for engineering smoke training and is therefore marked `development_exposed` for DBP1/DBP3. Folds 2-4 cover the previously unseen five proteins DBP5, DBP35, DBP48, DBP6, and DBP9. No model architecture, target definition, split, hyperparameter, or pair protocol was changed after the smoke result.

All metrics are calculated per protein on 8,192 canonical reverse-complement classes. Seed aggregation occurs within each protein/model before macro summarization. These results are not a row-level significance analysis.

## Primary all-7

| DBP | M0 | M1 | M1c | M2 | M3 | Delta M3-M1 | Delta M3-M1c | Delta M3-M2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DBP1 | 0.2223 | 0.3893 | 0.3691 | 0.0333 | 0.0500 | -0.3393 | -0.3191 | 0.0167 |
| DBP3 | 0.2236 | 0.0915 | 0.1346 | 0.1479 | 0.1119 | 0.0203 | -0.0227 | -0.0361 |
| DBP35 | 0.1250 | -0.0357 | -0.0130 | 0.0806 | 0.0477 | 0.0834 | 0.0606 | -0.0330 |
| DBP48 | -0.1585 | -0.1009 | -0.1554 | -0.1385 | -0.1951 | -0.0942 | -0.0397 | -0.0566 |
| DBP5 | 0.0324 | 0.0558 | 0.0605 | -0.0293 | 0.0682 | 0.0124 | 0.0077 | 0.0975 |
| DBP6 | 0.1967 | -0.0162 | 0.0052 | 0.0925 | 0.0724 | 0.0886 | 0.0673 | -0.0201 |
| DBP9 | 0.0962 | 0.2012 | 0.1842 | 0.0388 | 0.0715 | -0.1297 | -0.1127 | 0.0327 |

## Previously unseen five

| DBP | M0 | M1 | M1c | M2 | M3 | Delta M3-M1 | Delta M3-M1c | Delta M3-M2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DBP35 | 0.1250 | -0.0357 | -0.0130 | 0.0806 | 0.0477 | 0.0834 | 0.0606 | -0.0330 |
| DBP48 | -0.1585 | -0.1009 | -0.1554 | -0.1385 | -0.1951 | -0.0942 | -0.0397 | -0.0566 |
| DBP5 | 0.0324 | 0.0558 | 0.0605 | -0.0293 | 0.0682 | 0.0124 | 0.0077 | 0.0975 |
| DBP6 | 0.1967 | -0.0162 | 0.0052 | 0.0925 | 0.0724 | 0.0886 | 0.0673 | -0.0201 |
| DBP9 | 0.0962 | 0.2012 | 0.1842 | 0.0388 | 0.0715 | -0.1297 | -0.1127 | 0.0327 |

## Primary macro summary

| model | all7_macro_median | unseen5_macro_median | proteins_evaluated | unseen5_proteins_evaluated | seeds |
| --- | --- | --- | --- | --- | --- |
| M0 | 0.1250 | 0.0962 | 7 | 5 | 17|29|43 |
| M1 | 0.0558 | -0.0162 | 7 | 5 | 17|29|43 |
| M1c | 0.0605 | 0.0052 | 7 | 5 | 17|29|43 |
| M2 | 0.0388 | 0.0388 | 7 | 5 | 17|29|43 |
| M3 | 0.0682 | 0.0682 | 7 | 5 | 17|29|43 |

M3 improved over M1 on 4/7 proteins, over M1c on 3/7, and over M2 on 3/7. For the previously unseen five, the corresponding counts are 3/5, 3/5, and 2/5.

## Strict component sensitivity

This is the assay-informed conservative sensitivity split. It controls the combined protein-cluster/target-group/motif leakage components and is not the primary deployment estimate.

| DBP | M0 | M1 | M1c | M2 | M3 | Delta M3-M1 | Delta M3-M1c | Delta M3-M2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DBP1 | 0.2223 | 0.3893 | 0.3691 | 0.0333 | 0.0500 | -0.3393 | -0.3191 | 0.0167 |
| DBP3 | 0.2236 | 0.0915 | 0.1346 | 0.1479 | 0.1119 | 0.0203 | -0.0227 | -0.0361 |
| DBP35 | 0.1686 | -0.1159 | -0.1007 | 0.1216 | -0.1319 | -0.0160 | -0.0312 | -0.2535 |
| DBP48 | -0.1585 | -0.1009 | -0.1554 | -0.1385 | -0.1951 | -0.0942 | -0.0397 | -0.0566 |
| DBP5 | 0.0542 | -0.1487 | -0.1418 | 0.0133 | -0.0903 | 0.0584 | 0.0515 | -0.1035 |
| DBP6 | 0.2649 | 0.1478 | 0.1729 | 0.2377 | 0.1177 | -0.0300 | -0.0551 | -0.1200 |
| DBP9 | 0.0166 | -0.0134 | -0.0534 | -0.0203 | -0.0029 | 0.0105 | 0.0505 | 0.0174 |

| model | all_evaluated_macro_median | all_evaluated_macro_mean | all_evaluated_macro_sd | proteins_evaluated | components_evaluated | seeds |
| --- | --- | --- | --- | --- | --- | --- |
| M0 | 0.1686 | 0.1131 | 0.1510 | 7 | 3 | 17|29|43 |
| M1 | -0.0134 | 0.0357 | 0.1910 | 7 | 3 | 17|29|43 |
| M1c | -0.0534 | 0.0322 | 0.1976 | 7 | 3 | 17|29|43 |
| M2 | 0.0333 | 0.0564 | 0.1237 | 7 | 3 | 17|29|43 |
| M3 | -0.0029 | -0.0201 | 0.1223 | 7 | 3 | 17|29|43 |

Strict M3 median deltas were -0.0160 versus M1, -0.0312 versus M1c, and -0.0566 versus M2. The direction is not uniformly preserved across proteins.

## Target-relative controls

These controls use only the independently sourced `primary_target`; they do not use PBM-derived motifs.

| DBP | Edit | Hamming | Kmer overlap |
| --- | --- | --- | --- |
| DBP1 | 0.0454 | 0.0241 | 0.1330 |
| DBP3 | 0.0389 | 0.0085 | 0.1026 |
| DBP35 | -0.0000 | -0.0277 | 0.1354 |
| DBP48 | 0.2829 | 0.2599 | 0.2073 |
| DBP5 | 0.0606 | 0.0330 | 0.2133 |
| DBP6 | 0.1315 | 0.0856 | 0.2028 |
| DBP9 | 0.0789 | 0.0853 | 0.1104 |

## Seed stability

The per-protein seed standard deviations are stored in `primary_per_protein_results.csv` and `strict_component_per_protein_results.csv`. The complete run health summary is:

| quantity | value |
| --- | --- |
| primary model-runs | 60 |
| strict model-runs | 45 |
| minimum prediction variance | 0.0026 |
| maximum NaN/Inf count | 0 |
| total training runtime (seconds) | 234.3024 |
| training seeds | 17|29|43 |

## Baseline context

Prior baselines are shown for context only. Coverage and training regime differ, so their macro values are not unconditional rankings.

| baseline | coverage | macro_median_spearman | split_training_regime | directly_matched | notes |
| --- | --- | --- | --- | --- | --- |
| sequence-only kmer3 RC-aware | 7/7 | 0.2321 | v0.3.1 designed uPBM; no training | True | Sequence-only proxy; not protein-conditioned. |
| SimpleProteinConditional | 7/7 designed; 9 natural test | 0.3616 | natural PBM train -> designed external; prior v0.4.1 | False | Prior natural-test macro median=0.301354; low-capacity baseline. |
| FrozenPLM | 7/7 designed | 0.1530 | natural PBM train -> designed external; prior v0.4.2 | False | Frozen ESM-2 baseline; prior result is not this matched v0.5 trainer. |
| DeepPBS | 2/7 designed | 0.1593 | official structure-aware diagnostic; prior v0.4.2 | False | Coverage-limited diagnostic; not comparable as a seven-protein estimate. |
| NA-MPNN diagnostic | 2/7 designed | -0.0407 | official structure-aware diagnostic; prior v0.4.2 | False | Limited diagnostic coverage and DBP48 overlap caveat. |
| Experimental replicate reference | 7/7 designed | 0.5914 | uPBM replicate agreement | False | Empirical reproducibility reference, not a strict theoretical ceiling. |
| v0.5 M0 | 7/7 designed | 0.1250 | v0.5 4-fold protein-cluster LOCO, seed-mean per protein | True | Matched primary result; see primary tables for development-exposed and unseen subsets. |
| v0.5 M1 | 7/7 designed | 0.0558 | v0.5 4-fold protein-cluster LOCO, seed-mean per protein | True | Matched primary result; see primary tables for development-exposed and unseen subsets. |
| v0.5 M1c | 7/7 designed | 0.0605 | v0.5 4-fold protein-cluster LOCO, seed-mean per protein | True | Matched primary result; see primary tables for development-exposed and unseen subsets. |
| v0.5 M2 | 7/7 designed | 0.0388 | v0.5 4-fold protein-cluster LOCO, seed-mean per protein | True | Matched primary result; see primary tables for development-exposed and unseen subsets. |
| v0.5 M3 | 7/7 designed | 0.0682 | v0.5 4-fold protein-cluster LOCO, seed-mean per protein | True | Matched primary result; see primary tables for development-exposed and unseen subsets. |

## Interpretation

The frozen primary result does not provide robust evidence that M3 improves target-conditioned ranking over both protein-only matched controls: its all-7 median is below M1c, and its improvement counts are not consistent across the primary or strict split. M3 also does not exceed M2 in macro median. The unseen-five subset shows a small positive median delta versus M1 and M1c but remains negative versus M2 and is based on only five proteins.

This report does not issue a final GO/NO-GO decision. Hard-case analysis and failure-resolution analysis are intentionally deferred to the next phase. The current result is a falsification-oriented benchmark of the frozen minimal model family, not evidence that a final proposed model should be implemented.

## Future hypotheses

- Analyze whether the limited M3 effect is concentrated in specific target/motif groups or sequence-distance regimes.
- Separate implementation/capacity effects from genuine target-dependent ranking effects in the planned hard-case analysis.
- Keep the strict component result as a leakage-sensitivity reference rather than replacing the primary LOCO result.
