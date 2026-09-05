# v0.5 Phase 5A Smoke Report

## Status

This is a local-interaction engineering smoke test, not a primary scientific
result. The registered fold was `protein_cluster_loco_fold_2`, with DBP5 and
DBP35 held out. It is now `development-exposed-for-local-model`; folds 3 and 4
remain untouched.

The fixed protocol was retained: seed 42, 512 within-protein pairs per
training protein, 18 epochs, the existing logistic pairwise ranking loss,
learning rate 0.01, weight decay 1e-4, and the existing RC-class evaluation.

## Protein sequence audit

All seven designed sequences were taken from the existing Nature
supplementary design table and treated as reported designed-construct
sequences. They are complete for the residue encoder, but the exact PBM
expression tag and construct-boundary metadata were not independently
resolved. Lengths are DBP1 58, DBP3 58, DBP35 63, DBP48 65, DBP5 63, DBP6 56,
and DBP9 56 residues.

## Representation

The cache uses frozen ESM-2 `esm2_t12_35M_UR50D`, layer 12, with 480
dimensions per residue. BOS and EOS representations were removed, leaving one
embedding row per input residue. The cache contains 419 rows across the seven
proteins. No ESM parameter was updated.

L1 and L2 use one scaled dot-product attention block over protein residues.
Candidate and target orientations are pooled by arithmetic mean, and target
windows are pooled by arithmetic mean across canonical RC windows. Attention
values are mechanism diagnostics only, not contact predictions.

## Smoke results

The old columns below are frozen v0.5 primary seed-mean context values. Local
columns are the single Phase 5A smoke run and must not be compared as a new
primary estimate.

| DBP | old M1c | old M2 | old M3 | L1 | L1c | L2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DBP35 | -0.0130 | 0.0806 | 0.0477 | -0.1919 | -0.1888 | -0.1872 |
| DBP5 | 0.0605 | -0.0293 | 0.0682 | -0.3357 | -0.3356 | -0.3380 |

All local models scored 8192/8192 RC units per test protein and had zero
NaN/Inf predictions. The local test macro median Spearman values were L1
`-0.2638`, L1c `-0.2622`, and L2 `-0.2626`.

## Shuffle and attention diagnostics

| Diagnostic | Result |
| --- | ---: |
| L1 shuffled-P prediction correlation, median | 1.0000 |
| L2 shuffled-P prediction correlation, median | 1.0000 |
| L2 shuffled-T prediction correlation, median | 1.0000 |

The registered unchanged threshold was 0.995. The local models therefore did
not show material candidate-ranking sensitivity to the shuffled conditioning
inputs on this fold. Attention entropy was high (approximately 0.999 to 1.000
after normalization) and candidate-to-candidate attention variation was small
(approximately 0.0002 to 0.0004). There was no numerical collapse, but the
near-uniform attention is weak evidence that the trained block learned useful
candidate-specific residue selection.

## Training health

| Model | First loss | Final loss | Train median Spearman | Test variance | NaN/Inf |
| --- | ---: | ---: | ---: | ---: | ---: |
| L1 | 0.6905 | 0.6932 | 0.0240 | 9.08e-06 | 0 |
| L1c | 0.6909 | 0.6931 | 0.0290 | 6.09e-07 | 0 |
| L2 | 0.6946 | 0.6931 | 0.0182 | 2.59e-05 | 0 |

The first ReLU implementation produced a constant-prediction training-health
failure and was not used as a scientific result. That run is retained in
`logs/v0_5_local/01_initial_relu_failure.txt`. Replacing local MLP activations
with smooth tanh was the only repair; all registered data, split, seed,
objective, sampling, optimizer, and epoch settings were unchanged.

## Interpretation

The implementation and RC contract pass, but this one-fold result does not
support the Phase 5A hypothesis. L1 did not improve over the frozen global
protein-candidate context, and L2 did not demonstrate target use. The shuffled
conditioning diagnostics remain effectively unchanged, while attention is
near-uniform. This is evidence against the current local block under the
frozen sparse ranking protocol, not proof that all residue-level or
protein-DNA interaction representations are invalid.

**Decision: STOP.**

Do not run Phase 5B full local-model evaluation from this smoke result. The
next bottleneck should be reviewed before changing architecture: the current
512-pair training signal and its ability to learn candidate-dependent
protein-residue interactions. No proposed model was implemented.

### Current exposure correction

The exposure statements in this Phase 5A report are historical snapshots.
After the complete v0.5 primary evaluation and later diagnostics, all seven
designed DBPs have been viewed during development. They must not be described
as untouched confirmatory proteins. The unchanged historical artifacts remain
valid for Phase 5A provenance; current status is in
`metadata/v0_5_transfer/exposure_manifest.csv`.
