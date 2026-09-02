# v0.4 Proposed Model Requirements

This document is requirements-only. v0.4 does not implement a final model.

## Requirements Driven by Current Failure Analysis

1. The future baseline arena needs assay-matched natural PBM/uPBM data before training a protein-conditioned model. Without that control, natural-to-designed performance drops would be confounded by assay shift.
2. The model must be evaluated as a per-protein RC-class ranking task. The best current sequence-only baseline is `sequence_kmer3` with macro median Spearman 0.232, below empirical replicate agreement.
3. The method should not require a public complex structure for every designed DBP unless structure generation/mapping is itself benchmarked. v0.4 found public structure support for only DBP35 and DBP48 in checked sources.
4. The model must report missing predictions as not evaluable. Missing DeepPBS or SimpleProteinConditionalBaseline scores must never be filled with zero.
5. The model should explicitly target sequence-vs-experiment disagreement candidates. NA-MPNN resolved 50/398 evaluable v0.3.1 disagreement candidates in this diagnostic setup.
6. DBP48 must be handled as a non-zero-shot diagnostic for NA-MPNN because 8TAC appears in NA-MPNN split files.

## Not Yet Justified

- A large target-anchored neural architecture.
- Cross-protein absolute E-score calibration.
- Claims that designed DBPs are OOD without assay-matched natural PBM/uPBM controls.
