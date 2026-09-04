# v0.4.2 Progress

Date: 2026-09-02

## Completed

- Audited the natural UniPROBE PBM benchmark for assay construct provenance.
- Confirmed the local provenance does not recover assay-aligned construct sequences for the 57 natural proteins.
- Kept `FULL_LENGTH_REFERENCE` as a sensitivity benchmark only.
- Recorded DeepPBS official Linux runtime status and preserved Docker/WSL constraints in provenance.
- Trained and evaluated a frozen ESM-2 `esm2_t12_35M_UR50D` protein-conditioned baseline.
- Added difficulty-factor, training-space similarity, motif-distance, disagreement-resolution, and common-hard diagnostics.

## Current numbers

- Natural proteins audited: 57
- Assay-aligned construct benchmark rows: 0
- FrozenPLM natural-test macro median Spearman: 0.316
- FrozenPLM designed-external macro median Spearman: 0.153
- SimplePC designed-external macro median Spearman: 0.362
- Best prior sequence-only designed macro median Spearman: 0.232
- Empirical replicate Spearman reference: about 0.591
- v0.3.1 disagreement candidates: 1,515
- Core common-hard high-experiment/low-all-core cases: 263

## Decision and limitations

The gate remains `WAIT - BENCHMARK STILL INCOMPLETE`. The main blockers are
the empty assay-aligned construct benchmark and the absence of a runnable
Linux-compatible DeepPBS execution. FrozenPLM did not improve designed ranking
over SimplePC or the best sequence-only proxy.
