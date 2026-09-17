# Benchmark Protocol

## Scope

v0.7 is a frozen evidence synthesis and failure analysis. It reads existing v0.3-v0.6 artifacts and does not train new architectures, search hyperparameters, increase pair counts, or inspect quantitative labels from an external confirmatory cohort.

## Primary axis

The primary performance unit is within-protein Spearman ranking over canonical reverse-complement DNA units. The v0.5/v0.6 matched benchmark uses the four-fold protein-cluster LOCO split, seeds 17/29/43, the frozen pairwise objective and seven designed proteins.

## Independent axes

Performance reports median held-out Spearman, per-protein values and seed stability. Conditionality separately reports protein/target shuffle correlation and effect size (1 - correlation). A model with high ranking but shuffle correlation above 0.995 is classified high/low performance plus collapsed conditioning, not successful specificity learning.

## Exposure

All seven GSE237017 designed DBPs are development-exposed. They are valid for debugging and benchmark diagnosis only. No result in this directory is independent external validation. The master table preserves protocol, exposure and comparability fields; unmatched natural-to-designed SimplePC values remain context, not bridge evidence.

## Reproducibility

Historical v0.5 files are protected by their frozen manifest. v0.7 output is additive under `results/v0_7_benchmark/`, with machine-readable source paths and status fields.
