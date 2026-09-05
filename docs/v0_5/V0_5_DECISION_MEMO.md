# v0.5 Decision Memo

## Gate

**MODIFY**

The decision is based on the frozen primary evaluation plus the predeclared
v0.4.2 hard-case resolution rule. The full candidate-level replay matched the
primary seed metrics within `9.714e-17`.

## Evidence

- M3 joint capacity-controlled/target-only unique wins: **99**
  candidates, represented proteins: `DBP1|DBP3|DBP35|DBP48|DBP5|DBP6`.
- All five current models unresolved on: **871** hard-case candidates.
- M3 hard-case resolution: `236/1,515`
  (`15.6%`).
- Median M3 prediction correlation after protein shuffle:
  `1.0000`.
- Median M3 prediction correlation after target shuffle:
  `0.9999`.
- Median residualized M3 association after controlling TargetKmerOverlap:
  `0.1003`.

## One primary next hypothesis

**The frozen global ESM protein representation and low-capacity conditioning
head are the primary bottleneck to test before changing the target-conditioned
concept.**

Rationale: replacing P with a deterministic different designed-protein
embedding changed the frozen M3 prediction by a median correlation of
`1.0000`, while replacing T changed it by
`0.9999`. The protein-aware heads also show a
training/test gap in `train_vs_test_performance.csv`, and M3 does not beat the
capacity-matched M1c or target-only M2 in the frozen primary result. Pair
coverage is a separate limitation: 512 pairs expose a median
`11.7%` of candidates and all ten rank
deciles, so it is recorded but is not selected as the primary next hypothesis.

Minimal falsification experiment: keep the target manifest, primary LOCO split,
RC semantics, training objective, and evaluation fixed; replace the single
global mean-pooled protein vector with a pre-registered residue-level or
protein-DNA local representation, without using designed test outcomes for
representation selection. Compare it against the frozen global-embedding M3
on the same folds and seeds.

Falsifier: if a protein representation with demonstrably non-constant
shuffled-P sensitivity still fails to improve M3 relative to M1c and M2 on
previously unseen proteins, and does not produce a stable capacity-controlled
hard-case advantage, then representation alone is not the primary bottleneck.

## Not concluded

The analysis does not prove that target conditioning is biologically invalid.
It shows that the current minimal FiLM model, frozen ESM global embedding,
ranking objective, and sparse pair protocol do not yet provide robust evidence
for a target-conditioned advantage.
