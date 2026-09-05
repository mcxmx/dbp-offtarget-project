# Phase 6A Dense Pair Supervision Falsification Pilot

## Status

This is a one-fold falsification pilot, not a primary generalization result.
The registered smoke fold is `protein_cluster_loco_fold_3`, with DBP48 held
out and DBP1, DBP3, DBP5, DBP35, DBP6, and DBP9 used for training. The fold
was the first fold consisting only of the still-untouched proteins at the
start of Phase 6A. DBP48 is therefore now development-exposed for dense
supervision; DBP6 and DBP9 remain untouched.

The model family, labels, RC-class units, optimizer, learning rate, weight
decay, epochs, seed, pairwise logistic loss, and DNA/protein/target
representations were unchanged. The only registered experimental variable was
pair-supervision density and coverage.

## Registered protocols and observed coverage

| Protocol | Pairs/protein | Median unique candidates | Median coverage | Rank deciles |
| --- | ---: | ---: | ---: | ---: |
| S512 | 512 | 964 | 11.77% | 10/10 |
| D4096 | 4,096 | 6,278 | 76.64% | 10/10 |
| D16384 | 16,384 | 8,192 | 100.00% | 10/10 |

The dense sampler maintained the registered 40%/35%/25% easy, medium, and
hard strata and had zero duplicate unordered pairs in the audit. S512 used
the existing frozen sampler. D4096 and D16384 used the deterministic
coverage-aware sampler defined before test evaluation.

## Smoke training results

The test protein is DBP48, so each test rho below is a single-protein
diagnostic. Values are not macro estimates.

| Model | Protocol | Train rho | DBP48 test rho | First loss | Final loss |
| --- | --- | ---: | ---: | ---: | ---: |
| M0 | S512 | 0.3031 | -0.2030 | 0.6921 | 0.6273 |
| M1c | S512 | 0.3116 | -0.1643 | 0.6964 | 0.6165 |
| M2 | S512 | 0.2627 | -0.1833 | 0.6943 | 0.6344 |
| M3 | S512 | 0.3058 | -0.2559 | 0.6966 | 0.6192 |
| M0 | D4096 | 0.3330 | -0.2734 | 0.6832 | 0.6607 |
| M1c | D4096 | 0.3194 | -0.1819 | 0.6849 | 0.6623 |
| M2 | D4096 | 0.3021 | -0.3287 | 0.6822 | 0.6634 |
| M3 | D4096 | 0.2889 | -0.2931 | 0.6805 | 0.6704 |
| M0 | D16384 | 0.2722 | -0.2245 | 0.6607 | 0.6468 |
| M1c | D16384 | 0.2215 | -0.1180 | 0.6648 | 0.6470 |
| M2 | D16384 | 0.2950 | -0.2533 | 0.6662 | 0.6478 |
| M3 | D16384 | 0.2742 | -0.2690 | 0.6625 | 0.6519 |

All 12 runs completed without NaN/Inf values or collapsed test prediction
variance. Total recorded runtime was approximately 634 seconds on the local
CPU environment. D4096 and D16384 increased optimizer steps from 432 to 3,456
and 13,824, respectively, for the six training proteins.

The S512 DBP48 values reproduce the frozen v0.5 fold-3/seed-17 results:
M0 `-0.2030`, M1c `-0.1643`, M2 `-0.1833`, and M3 `-0.2559`. During
pre-validation, an extra `torch.set_num_threads(2)` caused a non-identical
CPU replay. That setting was removed before accepting these results; the
replay is now enforced by a test.

## Supervision effect

Relative to S512, the DBP48 test-rho changes were:

| Model | D4096 - S512 | D16384 - S512 |
| --- | ---: | ---: |
| M0 | -0.0704 | -0.0215 |
| M1c | -0.0176 | +0.0463 |
| M2 | -0.1454 | -0.0700 |
| M3 | -0.0372 | -0.0131 |

Training ranking did not improve consistently. D4096 improved training rho
for M0 and M2 but not M1c or M3; D16384 reduced training rho for M0, M1c,
and M3 relative to S512. The held-out changes were mixed or negative, with
only M1c at D16384 showing a positive change, while remaining below zero.

## Conditioning-use diagnostic

The shuffle tests were inference-only and did not retrain any model. A
deterministic source protein/target from DBP5 was substituted into the DBP48
test input.

| Model | Protocol | Shuffle | Prediction correlation | Mean absolute score change |
| --- | --- | --- | ---: | ---: |
| M1c | S512 | protein | 0.8321 | 0.5091 |
| M1c | D4096 | protein | 1.0000 | 0.0300 |
| M1c | D16384 | protein | 1.0000 | 0.0315 |
| M2 | S512 | target | 0.8951 | 0.2387 |
| M2 | D4096 | target | 0.9962 | 0.0254 |
| M2 | D16384 | target | 0.9973 | 0.0388 |
| M3 | S512 | protein | 1.0000 | 0.0033 |
| M3 | S512 | target | 0.9681 | 0.1642 |
| M3 | D4096 | protein | 1.0000 | 0.0004 |
| M3 | D4096 | target | 0.9814 | 0.0415 |
| M3 | D16384 | protein | 1.0000 | 0.0033 |
| M3 | D16384 | target | 0.9810 | 0.1099 |

Dense supervision did not produce a systematic increase in conditioning
sensitivity. For M1c, M2, and M3 the shuffled-input prediction correlations
were generally closer to one at the dense levels than at S512. This is not
evidence that the models learned more useful protein/target dependence.

## Pre-registered decision

### H2: NOT SUPPORTED

The registered hypothesis predicted that denser coverage would improve
training ranking, preserve healthy predictions, increase conditioning-use
sensitivity, and preferably improve held-out ranking for the same model.
Healthy numerical execution and candidate coverage were achieved, but the
learning, conditioning-use, and held-out-ranking parts were not consistent.

**Interpretation: STOP.** Do not run the full dense four-fold experiment,
change to a listwise objective, increase pair counts again, or transfer dense
supervision to L1/L2 based on this pilot. Sparse pair coverage may still be a
limitation, but this experiment does not establish it as the explanation for
the cross-protein failure.

### Current exposure correction

This report preserves the exposure state known during Phase 6A. It is a
historical snapshot, not the current confirmatory status. The complete v0.5
primary evaluation subsequently exposed all seven designed DBPs, including
DBP6 and DBP9. No GSE237017 designed protein is currently untouched; see
`metadata/v0_5_transfer/exposure_manifest.csv`.

The next possible hypothesis is outside this phase: the remaining bottleneck
may be cross-protein generalization or the frozen global protein
representation, rather than the number of sampled ranking pairs.
