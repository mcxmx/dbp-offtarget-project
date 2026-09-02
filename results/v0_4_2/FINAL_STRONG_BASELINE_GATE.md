# v0.4.2 Final Strong Baseline Gate

Date: 2026-09-02

## Decision

WAIT - BENCHMARK STILL INCOMPLETE

## Evidence

1. Natural PBM construct audit recovered no assay-aligned construct sequences from the current local provenance for the 57 UniPROBE proteins.
2. FrozenPLMProteinConditionalBaseline improved the natural-test macro median Spearman to 0.316.
3. FrozenPLM designed-external macro median Spearman is 0.153.
4. The prior SimplePC designed-external macro median Spearman was 0.362.
5. The best prior designed sequence-only baseline remains 0.232.
6. Designed uPBM empirical replicate Spearman reference remains about 0.591.
7. DeepPBS official Linux example was not run on this host because Docker/Podman/installed WSL are unavailable.
8. NA-MPNN remains diagnostic only; the prior v0.4 result covered 2/7 proteins.
9. The 1,515 sequence-vs-experiment disagreement candidates from v0.3.1 are unchanged as the disagreement reference set.

## Interpretation

FrozenPLM is a real, frozen protein-conditioned baseline and is stronger than the sequence-only proxy family, but it is still far below empirical replicate agreement on designed uPBM. Because the assay-aligned natural construct benchmark remains empty and DeepPBS is still not runnable on this host, the project should not yet transition to the final proposed model implementation from this repository state alone.
