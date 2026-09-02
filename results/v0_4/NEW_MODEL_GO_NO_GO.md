# v0.4 New Model Go / No-Go

Decision: CONDITIONAL GO

## Rationale

The fixed v0.3.1 benchmark is usable for baseline arena work, and the sequence-only gap is clear: best sequence-only macro median Spearman is 0.232, versus empirical replicate agreement reference 0.591.

However, this is not a STRONG GO yet. DeepPBS was not fairly runnable in the current environment, and NA-MPNN was evaluable for only 2/7 designed DBPs. Five designed DBPs have no public structure/model found in the checked official sources. DBP48/8TAC also has a high overlap risk because `8tac` appears in NA-MPNN split files.

## Required Before Strong Claims

- Run DeepPBS in a supported Linux/container environment or explicitly document it as not comparable.
- Expand structure availability or use structure-free baselines so all seven designed DBPs can be evaluated.
- Add assay-matched natural PBM/uPBM training data before training the Tier 1 protein-conditioned baseline.
- Keep DBP48/8TAC separate from zero-shot claims because of the detected NA-MPNN split overlap.

## Gate Interpretation

CONDITIONAL GO means v0.4 should continue into stronger baseline work and a careful Tier 1 assay-matched training setup. It does not authorize claims that existing strong baselines have systematically failed across the full designed-DBP benchmark.

## Key Constraints

- No new final model should be proposed until baseline coverage is improved.
- Any future comparison must remain per-protein and RC-class grouped.
- uPBM E-score remains an experimental specificity ranking signal, not affinity or in vivo binding.
