# v1.3.1 Prospective Holdout Evaluation

Prediction freeze commit: `ff68394`. Labels were joined only after that commit. The prediction artifact's cached SHA256 is `002c8825679a4a918bf7eac04c06eb390381061a52ae4f5197f5416cbb8dacdf`; Stage 2 checked only its cached path, size, and mtime and did not re-hash it.

## Per-protein decomposition

| protein | global rho | position rho | identity residual rho | identity pairwise accuracy | identity permutation p |
|---|---:|---:|---:|---:|---:|
| DBP023 | 0.308 | 0.560 | 0.006 | 0.548 | 0.350 |
| DBP056 | 0.085 | 0.455 | 0.001 | 0.452 | 0.743 |
| DBP062 | -0.030 | 0.705 | -0.038 | 0.571 | 0.237 |

## Three-protein summary

- Global median Spearman: 0.085.
- Position-sensitivity median Spearman: 0.560.
- Within-position identity median Spearman: 0.001.
- Within-position pairwise median accuracy: 0.548 (chance 0.5).
- The protein-level sample size is three. Bootstrap intervals are descriptive and no mutation-row pseudo-replication is used for biological inference.

## Interpretation gate

The prospective direction is consistent with Scenario D: position performance exceeds residual identity performance and identity pairwise accuracy remains near chance. This supports the failure-mode decomposition, but n=3 is too small for a broad significance claim.

Competition replicate identity reproducibility remains unavailable because the official workbook contains published means only. Therefore this holdout result alone cannot distinguish model-specific identity failure from a weakly identifiable assay endpoint, and v1.4 GNN development is not yet justified.

The development/holdout comparison is rendered from `results/v1_3/figure2_revised_data.tsv` as `figures/v1_3/Figure2_revised_development_holdout.png` and PDF. Development and locked-holdout proteins are visually separated.
