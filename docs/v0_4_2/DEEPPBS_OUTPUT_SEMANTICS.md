# DeepPBS Output Semantics

Audit date: 2026-09-04

This document records the source-level interpretation used by the project
adapter. The checked upstream checkout is `external/deeppbs/DeepPBS` at commit
`8bfb211dd67f02877841f6f33aa493ddf7daedf9`.

## Base order and input sequence

In `deeppbs/dna_encodings.py`, function `seqToOneHot`, the mapping is
`A -> 0`, `C -> 1`, `G -> 2`, and `T -> 3`. The inverse mapping in
`oneHotToSeq` uses the same `A/C/G/T` order. The same file's `rcSeq`
implements the standard reverse complement.

## Meaning of `P`

In `run/predict.py`, the upstream code applies `torch.softmax` to each
ensemble model output, sums and averages the ensemble outputs, then splits the
result into two strand halves. The first half and the reversed second half are
averaged:

`(output[:L] + flip(output[L:])) / 2`

The saved `P` array is therefore a position-wise, strand-averaged
post-softmax base-probability profile produced by the DeepPBS ensemble. It is
not a Kd, affinity, probability of binding, or calibrated specificity score.

The saved `Seq` array is the hard input sequence tensor converted from the
upstream batch. It is not an independent prediction.

## Orientation and RC handling

The upstream model explicitly constructs a reverse-complement DNA representation
in `run/models/model_v2.py::forward`, and the prediction script combines the
two strand halves. The project adapter additionally applies the benchmark's
frozen RC-class rule: a canonical 7-mer is scored in both orientations and the
larger PWM-derived log-probability is retained. This produces one score per
canonical reverse-complement class and does not duplicate those classes as
independent experimental units.

## Project-side 7-mer scorer

For a candidate 7-mer `x` and every contiguous seven-position window `j` in the
DeepPBS output profile, the primary score is:

`sum_i log(P[j+i, base(x_i)] + 1e-9)`

The maximum over valid windows is used. The maximum of the forward candidate
and its reverse complement is then used for the canonical RC class. This rule
was fixed before reading the PBM correlations and is implemented in
`src/deeppbs_adapter.py`.

## Provenance boundary

The model inference and preprocessing were performed by the unchanged
upstream DeepPBS checkout in an Ubuntu VM. The project only performs input
preparation for the documented DBP48 helix-only repair, NPZ validation,
canonical RC enumeration, and fixed score mapping. All output files and
checksums are recorded in
`results/v0_4_2/tables/deeppbs_run_manifest_completed_v0_4_2.csv`.
