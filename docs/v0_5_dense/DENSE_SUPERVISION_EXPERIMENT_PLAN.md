# v0.5 Phase 6A: Dense Pair Supervision Falsification Pilot

## Hypothesis

The frozen v0.5 P-conditioned models may underuse the dense PBM landscape
because each training protein contributes only 512 sampled ranking pairs while
having 8192 canonical reverse-complement candidate units. This pilot changes
only supervision density and candidate coverage.

## Frozen controls

The intended target, PBM labels, canonical RC units, protein-cluster LOCO
split, global M0/M1c/M2/M3 architectures, protein embeddings, optimizer,
learning rate, weight decay, activation, dimensions, loss, seed, and 18
epochs remain frozen. L1/L2 are excluded because their Phase 5A smoke exposed
optimization/attention collapse.

The ranking loss remains the original within-protein logistic pairwise loss.
Listwise objectives are explicitly out of scope so that supervision density is
the only experimental variable.

## Registered protocols

| Protocol | Pairs per protein | Sampler |
| --- | ---: | --- |
| S512 | 512 | existing frozen deterministic sampler |
| D4096 | 4096 | deterministic coverage-aware sampler |
| D16384 | 16384 | deterministic coverage-aware sampler |

The dense sampler uses only training-protein PBM ranks. It does not use
genome candidates, target similarity, model errors, or held-out labels. It
keeps the 40% easy, 35% medium, 25% hard rank-difference strata and prevents
duplicate unordered pairs.

The experiment intentionally allows optimizer steps and runtime to increase
with supervision density. Those costs are recorded rather than normalized
away because the hypothesis concerns underuse of supervision.

## Smoke fold and exposure

The first fold in manifest order whose test proteins were all still untouched
at the start of Phase 6A is `protein_cluster_loco_fold_3`, holding out DBP48.
The registered seed is 17. DBP48 becomes
`development_exposed_for_dense_supervision` after this run; DBP6 and DBP9
remain untouched.

This is a one-fold falsification pilot, not a primary generalization estimate.
No protocol will be selected after inspecting test results. All three
protocols will be reported.

## Pre-registered interpretation

- **H2 supported:** dense levels consistently improve training ranking,
  increase healthy prediction/conditioning sensitivity, and preferably
  improve held-out ranking for the same model.
- **Partially supported:** training improves but held-out ranking does not;
  sparse supervision explains undertraining but not cross-protein
  generalization.
- **H2 not supported:** dense supervision does not materially improve training
  or conditioning use.

No full four-fold dense experiment, listwise loss, local-model dense run, or
genome mining follows automatically from this pilot.
