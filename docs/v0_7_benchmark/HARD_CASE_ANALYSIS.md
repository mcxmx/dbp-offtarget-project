# Hard Case Analysis

The reference is the pre-registered v0.3.1 sequence-vs-experiment disagreement set (1,515 candidates). No thresholds were selected after inspecting v0.5/v0.6 predictions. These are ranking discrepancies, not biological binding-failure labels. The `reference_population`, `denominator_n`, `fraction_of_reference_1515` and `fraction_of_denominator` columns distinguish registered 1,515-case subsets, the 871 all-model-failure subset, and inference-only shuffle rows.

| category | count | fraction_of_reference_1515 | fraction_of_denominator | denominator_n | n_proteins | proteins | experimental_score_median | target_hamming_median | target_kmer_overlap_median |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sequence_only_rescued | 164 | 0.1083 | 0.1083 | 1515 | 6 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6 | 0.2513 | 0.4286 | 0.1429 |
| conditional_m1c_rescued | 154 | 0.1017 | 0.1017 | 1515 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2435 | 0.5714 | 0.1429 |
| conditional_m2_rescued | 118 | 0.0779 | 0.0779 | 1515 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2286 | 0.4286 | 0.1250 |
| joint_conditional_rescued | 99 | 0.0653 | 0.0653 | 1515 | 6 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6 | 0.2257 | 0.5714 | 0.1250 |
| all_model_failure | 871 | 0.5749 | 0.5749 | 1515 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2464 | 0.5714 | 0.1667 |
| protein_sensitive_cases | 1 | NA | 0.0476 | 21 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | NA | NA | NA |
| target_sensitive_cases | 10 | NA | 0.4762 | 21 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | NA | NA | NA |
| motif_near_target | 144 | 0.0950 | 0.1653 | 871 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2800 | 0.7143 | 0.6000 |
| highly_dissimilar_off_target | 50 | 0.0330 | 0.0574 | 871 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2514 | 0.8571 | 0.7750 |
| high_experimental_low_predicted | 481 | 0.3175 | 0.5522 | 871 | 7 | DBP1|DBP3|DBP35|DBP48|DBP5|DBP6|DBP9 | 0.2700 | 0.5714 | 0.2000 |
| low_experimental_high_predicted | 0 | 0.0000 | 0.0000 | 871 | 0 |  | NA | NA | NA |

## Findings

The strongest reproducible failure class is all-model failure: 871/1,515 (57.5%) candidates remain unresolved by M0-M3. Sequence-only rescue and conditional-model rescue subsets exist, but no stable joint conditional advantage was established (the pre-registered M3 joint-control subset contains 99 candidates). The motif/distance/score strata are calculated within the 871 all-model-failure subset; their denominator is explicitly recorded and they are descriptive only. Existing rows show broad representation across all seven DBPs rather than one protein-specific failure. GC, distance, overlap and score summaries are descriptive only; no small-sample causal interpretation is made. Protein/target sensitivity rows are model-level shuffle diagnostics, not candidate-level biological labels.
