# v0.5 Model Decision Memo

## Current recommendation

**WAIT FOR STRUCTURE-COVERAGE CLOSURE BEFORE CLAIMING A FULL STRONG-BASELINE
COMPARISON; DO NOT START THE PROPOSED MODEL YET.**

DeepPBS now runs through the official Linux workflow and supplies a valid
structure-aware diagnostic for DBP35 and DBP48. However, only **2/7**
designed proteins are evaluable. The observed median Spearman (0.159)
cannot be compared as a complete seven-protein result with the
0.232 sequence-only, 0.362 SimplePC, or
0.591 replicate references.

## What is supported

- DeepPBS preprocessing and inference are reproducible in the Ubuntu VM.
- The fixed PWM-to-7-mer mapping and RC-class scorer are implemented.
- DBP35 and DBP48 can be included as a transparent structure-aware diagnostic.
- Five proteins remain structurally unevaluable without fabricating inputs.

## What is not supported

- No claim that DeepPBS is the strongest baseline.
- No clean zero-shot claim because homolog-level training overlap is unresolved.
- No claim that a low result is biological rather than structure/mapping-related
  for DBP48 or theoretical-model-related for DBP35.

## v0.5 hypothesis to test after coverage closure

If future structure coverage confirms that DeepPBS, SimplePC, FrozenPLM, and
sequence-only baselines leave a stable common failure set, the minimal
failure-driven hypothesis is a target-anchored differential ranker:
compare `S(P,T)` and `S(P,D)` directly rather than scoring `D` in isolation.
The first falsification experiment should use a protein-held-out, RC-safe,
per-protein ranking benchmark and compare target-relative margins against the
best existing baseline with bootstrap confidence intervals. This is a memo
only; no proposed model is implemented in this commit.
