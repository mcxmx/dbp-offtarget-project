# v0.4 Failure Resolution Definition

This document is predeclared before inspecting v0.4 failure counts.

## Starting Set

Use v0.3.1 sequence-vs-experiment disagreement candidates: per protein, processed uPBM E-score consensus in the top 5% and RC-aware Hamming similarity to the paper motif at or below the per-protein median.

## Resolution Rule

For any protein-conditioned baseline with predictions on that protein, a v0.3.1 disagreement candidate is considered resolved when the baseline ranks it in the top 10% of that protein's predicted scores.

This asks whether a protein-conditioned method elevates experimentally high-scoring sequences that a sequence-only Hamming proxy did not prioritize. It does not imply the method correctly models physical binding mechanism.

## Not-Evaluable Cases

If a baseline has no prediction for a protein, candidates from that protein are counted as not evaluable for that baseline, not unresolved.
