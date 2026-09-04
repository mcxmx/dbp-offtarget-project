# v0.4.2 Disagreement Resolution Protocol

Date: 2026-09-02

This protocol is frozen before recomputing the v0.4.2 diagnostic tables.
The outputs describe ranking agreement with the experimental uPBM landscape.
They do not establish a physical binding mechanism or a binding affinity.

## Reference candidate set

The reference set is exactly the 1,515 sequence-vs-experiment disagreement
candidates defined in v0.3.1:

- processed uPBM E-score is at or above the within-protein 95th percentile;
- RC-aware motif Hamming similarity is at or below the within-protein 50th
  percentile threshold used in `all_disagreement_candidate_counts.csv`;
- the unit is one canonical reverse-complement 7-mer class.

The v0.3.1 thresholds and counts are reused as provenance. They are not
re-estimated after inspecting v0.4.2 predictions.

## Resolution rule

All prediction scores are oriented so that higher means a higher predicted
specificity ranking. For each protein and method, prediction percentiles are
computed across all available 8,192 canonical RC classes.

A reference disagreement candidate is:

- `resolved` if a prediction exists and its within-protein prediction
  percentile is at least 0.90;
- `unresolved` if a prediction exists and its percentile is below 0.90;
- `not_evaluable` if that method has no prediction for the candidate's protein.

Resolution rate is `n_resolved / n_evaluable`, not
`n_resolved / n_total_candidates`. The table reports all three counts.
Missing predictions are never replaced with zero.

The 0.90 threshold is a predeclared ranking triage rule. It asks whether a
method elevates an experimentally high-scoring sequence into its top decile;
it is not a calibrated probability threshold.

## Core common-hard set

The primary common-hard set is defined only using methods with complete,
structure-free predictions on all seven designed proteins:

- sequence-only `kmer3` proxy;
- `SimpleProteinConditionalBaseline`;
- `FrozenPLMProteinConditionalBaseline`.

A row is `high_experiment_low_all_core` when:

- experimental percentile is at least 0.95; and
- each core method's prediction percentile is at most 0.50.

The reciprocal `low_experiment_high_all_core` set is also reported for
diagnostic symmetry. DeepPBS is not included because its v0.4.2 prediction
table is empty, and NA-MPNN is not included because it covers only 2/7
proteins. These cases are therefore not claims that every available
protein-conditioned method failed.

## Method labels

The labels used in the output are:

- `sequence_kmer3`: RC-aware sequence similarity proxy;
- `SimpleProteinConditionalBaseline`: frozen v0.4.1 baseline;
- `FrozenPLMProteinConditionalBaseline`: frozen ESM-2 v0.4.2 baseline;
- `NA-MPNN diagnostic`: partial structural diagnostic only;
- `DeepPBS`: not evaluable in the current host state.

The word `failure` in filenames is shorthand for a ranking discrepancy. The
scientific prose uses `sequence-vs-experiment disagreement` or
`experimental specificity ranking failure case`, not biological binding
failure.
