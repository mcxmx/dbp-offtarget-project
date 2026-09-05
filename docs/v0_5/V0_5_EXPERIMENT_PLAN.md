# v0.5 Experiment Plan

## Scope

v0.5 defines the data semantics and split contract for a future
target-conditioned specificity model. This document does not train or
implement that model.

The intended task is:

```text
R(P, T, D)
```

where `P` is a protein, `T` is an independently reported intended/design
target, and `D` is a candidate DNA sequence. The primary target is never
selected from PBM top-scoring sequences, PBM E-scores, or a test-set optimum.

The project keeps three distinct DNA concepts:

1. `intended_design_target`: sequence reported for the design experiment.
2. `experimental_assay_reference`: DNA sequence used as a later assay/reference
   where the provenance distinguishes it from the design target.
3. `pbm_motif`: shorter motif used to interpret the PBM landscape. It is not a
   replacement for the intended target.

DBP48 is the explicit audit case: its original design target is
`CGCCCAAAGCCGCG`, its experimental/PBM reference is
`CGACACCTGACGCG`, and its PBM motif is `CTGACG`.

## Why an intended target must be part of the task

A separable score of the form

```text
S(P, T) - S(P, D)
```

does not solve the target-conditioned ranking problem when `T` is fixed for a
protein. The first term is a candidate-independent constant. For any two
candidates `D1` and `D2`:

```text
S(P,T) - S(P,D1) > S(P,T) - S(P,D2)
iff
S(P,D1) < S(P,D2)
```

Therefore subtracting `S(P,T)` cannot change the within-protein ranking of
candidate DNA sequences and cannot improve within-protein Spearman by itself.

The proposed model must be non-separable:

```text
R(P, T, D)
```

Here `T` must change how the model interprets `D`, for example through a
comparative representation, target-relative interaction, or another mechanism
that is tested rather than assumed.

## Split contract

The PBM-derived motif may be used for the conservative leakage-sensitivity
grouping recorded in the v0.5 audit. It must not be used as the primary target,
as a model input, as feature engineering, or for test-performance tuning.

### Primary split: leave-one-protein-cluster-out

The primary unit is a whole designed protein sequence cluster, using the
existing v0.3.1 cluster definitions. No protein may occur in both train and
test within a fold. The current four folds are:

- `protein_cluster_1`: DBP1, DBP3
- `protein_cluster_2`: DBP35, DBP5
- `protein_cluster_3`: DBP48
- `protein_cluster_4`: DBP6, DBP9

This is the primary protein-level stress test, but it is not fully target-group
independent because the TGCACA motif is shared across protein clusters 2 and 4.

### Strict sensitivity split: combined leakage-component-out

The leakage graph connects proteins sharing a protein cluster, original target
group, assay target group, or motif group. The current graph produces three
components:

- DBP1, DBP3
- DBP35, DBP5, DBP6, DBP9
- DBP48

The component-out split is the strict sensitivity analysis. It has only three
folds and therefore has limited statistical resolution, but it prevents the
known target/motif connections from crossing train and test.

### DNA unit rule

Future DNA-level evaluation uses `canonical_rc_equivalence_class` as the unit.
An oriented sequence and its reverse complement may not be placed in different
partitions. Randomly splitting the 57,344 protein-RC-class rows is prohibited.

## Planned evaluation

The future model should report per-protein and per-fold ranking metrics first,
then macro summaries with uncertainty. The designed DBPs should remain an
external/designed benchmark for models trained on natural specificity data
unless a clearly labeled within-designed diagnostic is being run.

No model, hyperparameter, epoch, or checkpoint may be selected using the
designed test folds.

## Development Exposure Note

The engineering smoke test used `protein_cluster_loco_fold_1`, whose held-out
proteins are DBP1 and DBP3. The architecture, target definition, split
contract, pair-sampling protocol, and hyperparameters were not changed in
response to those outputs, but this fold is nevertheless
`development-exposed` rather than untouched confirmatory evidence.

The previously unseen primary folds are folds 2-4, covering DBP5, DBP35,
DBP48, DBP6, and DBP9. The complete seven-protein result will report all
folds, while the primary confirmatory descriptive summary will also report
this previously unseen five-protein subset separately.

## Frozen Evaluation Seeds

Before evaluating folds 2-4, the multi-seed set is fixed to `17`, `29`, and
`43`. Seed aggregation is performed within each DBP and model first
(mean, standard deviation, minimum, and maximum), followed by macro summaries
over proteins. The seven proteins are not treated as 21 independent
biological replicates, and seed variation is not biological uncertainty.
