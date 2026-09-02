# v0.4.1 Final Model Development Gate

Date: 2026-09-02

## Decision

WAIT FOR STRONGER BASELINE

## Evidence

1. Natural PBM benchmark: 57 proteins and 1,875,072 protein-RC-class 8-mer units are now available for simple training.
2. Natural held-out SimplePC macro median Spearman: 0.301.
3. Designed external SimplePC macro median Spearman: 0.362.
4. Best prior sequence-only designed baseline macro median Spearman: 0.232.
5. Designed uPBM empirical replicate Spearman reference: about 0.591.
6. DeepPBS was not fairly run in v0.4.1 because this host lacks Docker/WSL runtime for the official Linux preprocessing workflow.
7. NA-MPNN remains diagnostic only: 2/7 designed proteins covered; DBP48 has known overlap risk.
8. SimplePC resolves 333/1,515 v0.3.1 sequence-vs-experiment disagreement candidates (22.0%), leaving most unresolved.

## Interpretation

SimplePC improves over the sequence-only baseline on designed DBPs, but it remains well below the empirical replicate reference and fails on DBP6/DBP48. Because DeepPBS has not yet been fairly reproduced, the project should not start the final proposed model as the next step. The immediate next baseline task is a Linux/Docker DeepPBS run on all structurally evaluable designed DBPs, followed by the same per-protein ranking evaluation.
