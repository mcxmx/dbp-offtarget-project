# Preregistered Validation Plan (Template; Labels Locked)

**Status:** FROZEN TEMPLATE, not executed. No external quantitative labels have been read.

1. **Primary hypothesis:** frozen methods have measurable transferable DNA ranking on an independent designed-DBP cohort under strict unseen-protein evaluation.
2. **Methods:** sequence k-mer, frozen v0.6 M0-M3, SimplePC/FrozenPLM context only, and any structure method with verified compatible inputs; no fine-tuning.
3. **Checkpoint/code:** hashes are recorded in `FROZEN_VALIDATION_MANIFEST.yaml`.
4. **Scoring:** canonical 7-mer within-protein ranking; fixed PWM/PPM log-probability conversion for structure methods.
5. **Register/strand:** predeclared biological register; reverse-complement classes are scored invariantly; ambiguous register means exclusion, never best-result alignment.
6. **Primary metric:** median per-protein Spearman.
7. **Secondary metrics:** per-protein Spearman, top-k/NDCG where pre-specified, hard-case recovery, and conditionality where meaningful.
8. **Exclusions:** missing construct/target, incompatible assay, ambiguous register, or incomplete scores according to predeclared thresholds.
9. **Missing data:** retain all eligible proteins in denominator; report coverage and failures.
10. **Failure handling:** preserve failed outputs and logs; no silent retries with changed settings.
11. **No post-hoc selection:** no model, checkpoint, register, metric, or protein removal after label access.
12. **Statistics:** protein is the independent unit; paired exact/permutation summaries and effect sizes; no candidate-level significance claims.
13. **Confirmatory vs exploratory:** the primary frozen pipeline is confirmatory; every later analysis is labeled POST-HOC / EXPLORATORY.

External validation remains blocked until a cohort satisfies the selection protocol.
