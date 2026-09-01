# Sequence-vs-Experiment Disagreement Definition

Audit date: 2026-09-01

The v0.3 file `sequence_similarity_failure_candidates.csv` contained 140 rows because the analysis selected at most 20 example rows per protein. That number is an example-table size, not the total number of discoveries.

## v0.3.1 Terminology

Use: `sequence-vs-experiment disagreement cases`.

Do not use: `model failures`.

Reason: no protein-conditioned model has been evaluated yet.

## v0.3.1 Criterion

For each protein separately, using 8192 reverse-complement equivalence classes:

- Experimental high-score condition: processed uPBM E-score consensus >= the per-protein 95% quantile.
- Sequence-similarity condition: RC-aware hamming similarity to the paper motif <= the per-protein 50% quantile.

This is a percentile-based descriptive rule. It is not an assay-defined binding threshold and should not be interpreted as an absolute off-target risk threshold.

## Outputs

- `results/v0_3_1/tables/all_disagreement_candidate_counts.csv`: total counts by protein.
- `results/v0_3_1/tables/top_disagreement_examples.csv`: at most 20 examples per protein, for inspection only.
