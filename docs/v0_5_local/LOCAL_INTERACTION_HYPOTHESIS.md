# v0.5 Phase 5A: Local Interaction Hypothesis

## Scope

This phase tests one representation hypothesis only:

> A single global mean-pooled protein embedding does not provide enough
> protein-specific information for cross-protein DNA specificity ranking.
> Candidate DNA should interact with frozen per-residue protein embeddings so
> different candidates can select different residue features.

This is a pilot, not a proposed final model. No ESM parameter is fine-tuned,
and no target definition, split, ranking loss, pair sampler, optimizer, seed,
or epoch count is changed from the frozen v0.5 protocol.

## Frozen design

| Item | Fixed choice |
| --- | --- |
| Protein encoder | frozen ESM-2 `esm2_t12_35M_UR50D`, layer 12, 480 dimensions |
| Residue tokens | BOS/EOS removed; one embedding per input residue |
| Local interaction | one scaled dot-product residue attention block |
| Local hidden dimension | 24 |
| Activation | smooth `tanh` in the small interaction and scoring MLPs |
| Candidate RC rule | mean of sequence and reverse-complement oriented scores |
| Target windows | all canonical RC 7-mer windows from the primary target |
| Target-window RC rule | mean over both orientations |
| Target-window pooling | arithmetic mean |
| Loss | frozen within-protein logistic pairwise ranking loss |
| Pair protocol | 512 pairs/protein, existing 40/35/25 difficulty strata |
| Optimizer | existing Adam, learning rate 0.01, weight decay 1e-4 |
| Epochs | 18 |
| Seed | 42 |

The capacity control `L1c` increases only the P,D head width to 64. Its
purpose is to keep the P,D-only local control in the same parameter range as
L2; it does not receive target features.

## Smoke fold registration

The first previously-unseen primary fold in the frozen manifest is
`protein_cluster_loco_fold_2`, holding out DBP5 and DBP35. It is registered as
the Phase 5A smoke fold before local-model training and will be labeled
`development-exposed-for-local-model` in all reports. DBP1/DBP3 were already
development-exposed by the earlier global-model smoke test. Folds 3 and 4
remain unrun in this phase.

The smoke seed is fixed to 42 because this is the existing v0.5 engineering
seed. This phase does not run a full primary evaluation or use test scores for
model selection.

## Pre-registered interpretation

- If L1 improves over the frozen global P,D control and shuffled protein
  materially changes predictions, this supports global representation as a
  bottleneck.
- If L1 is similar to the global control and shuffled protein is nearly
  unchanged, residue-level representation alone does not solve the issue.
- If L2 improves over L1/L1c and shuffled protein and target both change
  predictions, this is evidence to justify a Phase 5B full local-model
  evaluation.
- A low smoke Spearman is not a technical failure unless prediction variance
  collapses, NaN/Inf values occur, or the run violates the split/RC contract.

The material shuffle threshold is fixed before training: a median prediction
correlation greater than `0.995` is treated as effectively unchanged for this
pilot. Attention values are mechanism diagnostics only and are not interpreted
as biological contact maps.

The first implementation smoke exposed a constant-prediction failure in the
local ReLU head at the inherited learning rate. This was treated as an
implementation stability failure, not a scientific result; the local MLP
activation was changed to smooth `tanh` before the replacement smoke, with all
training and evaluation controls above unchanged.
