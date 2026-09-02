# v0.4.1 DeepPBS Task Mapping

DeepPBS predicts DNA specificity profiles/PWMs from a protein-DNA complex structure after Linux-oriented structural preprocessing. The designed-DBP benchmark target is a per-protein ranking of 7-mer reverse-complement classes by processed uPBM E-score.

Pre-registered mapping for a future runnable Linux execution:

1. Use only author-released or experimental protein-DNA structures listed in `metadata/v0_4_1/deeppbs_structure_manifest.csv`.
2. Run the official DeepPBS preprocessing and ensemble inference without modifying third-party model code.
3. Convert the predicted PWM/profile to candidate 7-mer scores by summing log probabilities across the aligned structure-defined DNA window.
4. Evaluate only per protein with Spearman, NDCG@1%, NDCG@5%, pairwise ranking accuracy, and top-1% experimental recovery.

The PWM-derived score must not be described as PBM affinity, Kd, or calibrated binding probability. If fewer than four designed DBPs have structures and runnable preprocessing, DeepPBS remains a partial diagnostic baseline.
