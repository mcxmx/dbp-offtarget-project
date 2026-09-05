# v0.5 Phase 2 Smoke Report

## Scope

This report covers implementation validation only. It uses the first
predefined primary `protein_cluster_loco` fold and is explicitly **NOT PRIMARY
SCIENTIFIC RESULT**. The complete four-fold evaluation was not run.

## Fixed fold

- Split: `protein_cluster_loco`
- Fold: `protein_cluster_loco_fold_1`
- Training proteins: DBP35, DBP48, DBP5, DBP6, DBP9
- Held-out proteins: DBP1, DBP3
- Units per protein: 8,192 canonical reverse-complement classes
- Training pairs: 2,560 total, 512 per training protein
- Seed: 42
- Pair strata: 40% easy, 35% medium, 25% hard

No held-out protein labels were used for pair selection, training, or
configuration selection.

## Model and parameter audit

| Model | Inputs | Trainable parameters |
|---|---|---:|
| M0_CandidateDNAOnly | D | 2,017 |
| M1_ProteinCandidate | P,D | 19,457 |
| M1c_ProteinCandidateCapacityMatched | P,D | 25,281 |
| M2_TargetCandidateOnly | T,D | 6,145 |
| M3_ProteinTargetCandidate | P,T,D | 23,649 |

M1 and M3 use the same external frozen ESM-2
`esm2_t12_35M_UR50D` mean-pooled 480-dimensional representation. The protein
embedding is not trainable. M3 changes the comparative target/candidate hidden
representation through protein-controlled FiLM modulation.

## Smoke ranking metrics

| DBP | M0 | M1 | M2 | M3 | M1c |
|---|---:|---:|---:|---:|---:|
| DBP1 | 0.2095 | 0.3148 | -0.1321 | 0.2494 | 0.3903 |
| DBP3 | 0.2367 | 0.1179 | 0.0785 | 0.1462 | 0.1895 |

These values demonstrate that the end-to-end training and evaluation path
produces finite, non-constant predictions. They must not be interpreted as a
primary model comparison because only one fold was run.

## Training health

All five models completed 18 epochs. Pairwise loss decreased for every model:

| Model | First loss | Last loss | Prediction variance | NaN/Inf |
|---|---:|---:|---:|---:|
| M0 | 0.7024 | 0.6450 | 0.2060 | 0 |
| M1 | 0.7036 | 0.5831 | 0.9833 | 0 |
| M2 | 0.7091 | 0.6551 | 0.5166 | 0 |
| M3 | 0.7101 | 0.6050 | 0.2310 | 0 |
| M1c | 0.7125 | 0.5828 | 0.4997 | 0 |

Total wall-clock training time for the five models was approximately 7.2
seconds in the local Python 3.13 environment.

## RC and non-separability checks

- Candidate sequence and its reverse complement receive identical scores for
  every model within numerical tolerance.
- Reversing the intended target orientation leaves M3 unchanged.
- Changing `T` changes M3 output for fixed `P,D`.
- Changing `P` changes M3 output for fixed `T,D`.
- A toy interaction test demonstrates that changing `T` can reverse candidate
  ranking, unlike adding a candidate-independent target constant.

## Target-blind controls

On the held-out DBP1/DBP3 smoke fold, Spearman values were:

| Control | DBP1 | DBP3 |
|---|---:|---:|
| TargetHamming | 0.0241 | 0.0085 |
| TargetEdit | 0.0454 | 0.0389 |
| TargetKmerOverlap | 0.1330 | 0.1026 |

These are protein-blind sequence controls, not binding-specificity models.

## Conclusion

The minimal v0.5 matched model family, RC contract, deterministic pairwise
trainer, parameter audit, and one legal smoke fold are operational. The
implementation is ready for review before any complete four-fold primary
experiment. No architecture or hyperparameter was selected using designed
held-out results, and no proposed-model claim is made here.
