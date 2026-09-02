# v0.4.2 Progress

Date: 2026-09-02

## Completed

- Audited the natural UniPROBE PBM benchmark for assay construct provenance.
- Confirmed the local provenance does not recover assay-aligned construct sequences for the 57 natural proteins.
- Kept `FULL_LENGTH_REFERENCE` as a sensitivity benchmark only.
- Recorded DeepPBS official Linux runtime status and preserved Docker/WSL constraints in provenance.
- Trained a frozen ESM-2 `esm2_t12_35M_UR50D` protein-conditioned baseline with a small ridge head.
- Evaluated the frozen baseline on natural held-out and designed external sets.
- Wrote the v0.4.2 summary gate and validation report.

## Current numbers

- Natural proteins audited: 57
- Assay-aligned construct benchmark rows: 0
- FrozenPLM natural-test macro median Spearman: 0.316
- FrozenPLM designed-external macro median Spearman: 0.153
- SimplePC designed-external macro median Spearman: 0.362
- Best prior sequence-only designed macro median Spearman: 0.232
- Empirical replicate Spearman reference: about 0.591

## Limitation

The construct-aware natural PBM benchmark remains empty, so the current protein-conditioned comparison is still anchored to full-length reference sequences rather than verified assay constructs.
