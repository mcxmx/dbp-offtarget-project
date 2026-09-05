# v0.5 Phase 2 Model Implementation

## Scope

This phase implements a matched model family and one engineering smoke test.
It does not run the complete four-fold primary evaluation and does not claim a
scientific performance result.

The models are:

| Model | Inputs | Role |
|---|---|---|
| `M0_CandidateDNAOnly` | `D` | DNA-only baseline |
| `M1_ProteinCandidate` | `P,D` | protein-conditioned candidate baseline |
| `M1c_ProteinCandidateCapacityMatched` | `P,D` | M1 capacity control |
| `M2_TargetCandidateOnly` | `T,D` | protein-blind target-relative control |
| `M3_ProteinTargetCandidate` | `P,T,D` | minimal non-separable target-conditioned model |

`M3` is a baseline/prototype implementation, not the final proposed method.

## Frozen representations

M1, M1c, and M3 use the same precomputed frozen ESM-2
`esm2_t12_35M_UR50D` mean-pooled representation already present in
`data/interim/v0_4_2/frozen_plm_embeddings_esm2_t12_35M_UR50D.parquet`.
The representation is 480-dimensional and is not fine-tuned.

Candidate DNA uses RC-symmetric one-hot 7-mer features. For every sequence,
the one-hot vectors for `D` and `RC(D)` are averaged before the learned
projection. Target sequences are represented by all unique canonical RC
7-mer windows and mean pooled. This makes target orientation arbitrary.

## M3 non-separability

M3 computes candidate/target comparative features:

```text
concat(d, t, d - t, d * t)
```

and applies a protein-controlled FiLM transformation before the scalar head.
Protein therefore changes the comparative hidden representation, rather than
being added as a candidate-independent offset at the end.

The algebraic test in `tests/test_v0_5_models.py` separately verifies that a
constant target term cannot reverse candidate ranking, while a toy
non-separable interaction can.

## Training contract

- Fixed seed: `42`.
- Training pairs are sampled only from the train proteins in
  `metadata/v0_5/v0_5_split_manifest.csv`.
- Pairs are within-protein and deterministically sampled into easy/medium/hard
  rank-difference bins.
- Primary objective: logistic pairwise ranking loss.
- No test protein contributes labels, pairs, validation feedback, or
  hyperparameter selection.
- Smoke fold: first deterministic primary fold,
  `protein_cluster_loco_fold_1`.
- Smoke output is explicitly `NOT PRIMARY SCIENTIFIC RESULT`.

## Target controls

`TargetHamming`, `TargetEdit`, and `TargetKmerOverlap` are protein-blind
sequence controls. They compare each candidate against the independently
reported `primary_target` windows, with both orientations considered. They do
not use the PBM-derived motif.
