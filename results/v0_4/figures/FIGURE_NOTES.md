# v0.4 Figure Notes

All v0.4 panels evaluate processed experimental uPBM E-score ranking, not binding affinity or off-target prediction.

- `fig_v0_4_1_baseline_performance_overview.png`: macro median per-protein Spearman. This shows broad ranking performance and evaluation coverage; it does not prove absolute calibration.
- `fig_v0_4_2_per_protein_spearman_heatmap.png`: per-DBP Spearman. Gray/NA cells are not-evaluable methods, not zero-valued failures.
- `fig_v0_4_3_prediction_vs_experimental_examples.png`: DBP35 and DBP48 examples for sequence-only and NA-MPNN structural PPM scores. These are diagnostic mappings to PBM 7-mer ranking.
- `fig_v0_4_4_performance_vs_replicate_reference.png`: compares evaluated baselines with empirical replicate agreement. Replicate agreement is an assay reproducibility reference, not a strict upper bound.
- `fig_v0_4_5_performance_by_motif_distance.png`: checks whether performance changes with RC-aware motif-distance regime.
- `fig_v0_4_6_failure_case_landscape.png`: separates total sequence-vs-experiment disagreement cases from examples resolved by NA-MPNN where predictions exist.
