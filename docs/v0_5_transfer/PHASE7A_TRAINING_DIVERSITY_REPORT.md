# v0.5 Phase 7A Training-Diversity Report

## Decision

**H3: NOT SUPPORTED**

The natural-transfer mainline is stopped. No new transfer model was trained,
no designed labels were added to natural training, and no target-conditioned
architecture was implemented.

## Exposure correction

The complete v0.5 primary evaluation viewed all seven designed proteins:
DBP1, DBP3, DBP5, DBP35, DBP48, DBP6, and DBP9. They are all
`development_exposed`. There is no independent confirmatory GSE237017 subset
remaining in this repository.

The historical smoke-phase fields are preserved unchanged. The current
exposure state is in `metadata/v0_5_transfer/exposure_manifest.csv`.

## Evidence boundary

The existing prior SimplePC values (`0.301` natural held-out and `0.362`
designed external) are retained as unmatched context only. The repository
does not contain a valid Phase 7A replay or a causal comparison of
protein-diversity versus supervision-quantity regimes.

Accordingly, this phase does not claim that natural training diversity was
tested and supported. It records that the H3 mainline was stopped before
bridge execution, as requested.

## What remains unknown

- Whether a same-architecture natural-only transfer model would outperform
  designed-only training under the frozen v0.5 contract.
- Whether any gain would come from protein diversity, additional observations,
  assay prior, or protocol differences.
- Whether exact experimental construct sequences would change the transfer
  result.

These are unanswered rather than filled with synthetic results.

## Next single falsifiable hypothesis

**H4: assay/task alignment, rather than protein diversity alone, is the
dominant transferable-learning bottleneck.**

The current natural benchmark uses processed 8-mer UniPROBE scores and mostly
reference full-length protein sequences, whereas the designed benchmark uses
processed 7-mer GSE237017 uPBM E-scores. A small, independently curated
natural PBM set with verified experimental constructs and matched k-mer/score
processing should be tested with the same frozen P-D architecture and
protein-level split contract.

Minimal falsification experiment: construct that assay-matched benchmark with
pre-registered train/validation/test proteins, train the already existing
bridge architecture without designed test labels, and evaluate the
developmental designed set with protein-shuffle diagnostics.

Falsifier: if the designed performance gap persists under matched assay
processing and verified construct sequences, while protein shuffling changes
the ranking, then assay/task alignment is not the dominant explanation.

This hypothesis is recorded only. It is not implemented in Phase 7A.

## Validation requirement

Any future selected model requires an independent designed-DBP dataset or a
prospective holdout that was not viewed during development.
