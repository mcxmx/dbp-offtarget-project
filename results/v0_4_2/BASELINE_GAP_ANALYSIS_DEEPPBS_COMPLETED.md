# Baseline Gap Analysis After DeepPBS Completion

Audit date: 2026-09-04

## Main result

DeepPBS is technically reproducible in the Ubuntu VM, but it is evaluable for
only **2/7** designed DBPs. Its diagnostic median Spearman is
**0.159271**, compared with the frozen sequence-only best baseline
**0.232082**, SimpleProteinConditional **0.361611**,
FrozenPLM **0.153020**, and empirical replicate agreement
**0.591439**.

These are not directly comparable as a complete model ranking because coverage,
structure requirements, and training/overlap caveats differ.

## Interpretation

- DeepPBS does not demonstrate a clear improvement over the existing
  sequence-only or SimplePC results on the two evaluable proteins.
- DBP35 is low (`0.040176`)
  and DBP48 is higher (`0.278366`),
  but the sample size is two and DBP35 uses a theoretical design model.
- The previous 1,515 disagreement count remains an all-protein count.
  DeepPBS has 398 eligible candidates, resolves 39, and leaves 359 unresolved.
- The completed four-method common-low set contains
  25 candidates.

## Boundary of claim

This closes the runtime and evaluation gap, but not the structure-coverage gap.
The result supports a structure-aware diagnostic baseline and a reproducible
failure analysis. It does not support a clean zero-shot or seven-protein
generalization claim.
