# V0.6 Results

## Status and scope

This is one frozen replay after representation and synthetic tests passed.
The replay uses the v0.5 protein-cluster LOCO split, seeds `17|29|43`, 512
within-protein pairs, the same logistic pairwise objective, 18 epochs, and no
hyperparameter search. The seven designed DBPs are all development-exposed;
these results are not independent external validation.

## Historical evidence boundary

The repair is interpreted together with the already frozen v0.5 diagnostics:
871/1,515 pre-registered hard cases were unresolved by every M0-M3 model and
M3 resolved only 236/1,515. The v0.5 M3 protein-shuffle and target-shuffle
medians were already 1.0000 and 0.9999. The dense-supervision H2 pilot did not
improve held-out ranking or conditioning sensitivity, and the local/residue
H1 pilot had near-uniform attention with protein and target shuffle medians at
1.0000. Phase 7A H3 natural-to-designed transfer was explicitly stopped as
`NOT_SUPPORTED`; its unmatched prior SimplePC numbers are not a bridge result.
These constraints are why v0.6 performs a representation audit and one fixed
replay rather than further pair scaling, attention stacking, or transfer
tuning.

## v0.5 versus v0.6

| model | v0.5 all-7 median | v0.6 all-7 median | delta |
| --- | ---: | ---: | ---: |
| M0 | 0.1250 | 0.2829 | 0.1580 |
| M1 | 0.0558 | 0.2016 | 0.1458 |
| M1c | 0.0605 | 0.1693 | 0.1088 |
| M2 | 0.0388 | 0.1342 | 0.0955 |
| M3 | 0.0682 | 0.1415 | 0.0733 |

v0.6 median Spearman is higher for M0/M1/M1c/M2/M3 respectively by the values
in `v05_vs_v06_comparison.csv`, but this is a developmental replay over an
already exposed cohort. It is evidence that the RC repair changes the learned
landscape, not evidence of confirmatory generalization.

## v0.6 per-protein Spearman (seed mean)

| DBP | M0 | M1 | M1c | M2 | M3 |
| --- | ---: | ---: | ---: | ---: | ---: |
| DBP1 | 0.3324 | 0.4137 | 0.4484 | -0.1227 | 0.0421 |
| DBP3 | 0.2829 | 0.2232 | 0.2251 | 0.1480 | 0.2099 |
| DBP5 | 0.2092 | 0.2490 | 0.1658 | 0.1271 | 0.1415 |
| DBP35 | 0.2841 | 0.1542 | 0.0871 | 0.2209 | 0.1792 |
| DBP48 | 0.0001 | -0.0160 | 0.0742 | -0.0051 | -0.1076 |
| DBP6 | 0.3187 | 0.1371 | 0.1693 | 0.1342 | 0.0893 |
| DBP9 | 0.1536 | 0.2016 | 0.2022 | 0.1684 | 0.1578 |

The sequence-only k-mer baseline remains a historical v0.3.1 comparator with
all-7 median Spearman `0.2321`; it is not retrained in v0.6. v0.6 medians are
M0 `0.2829`, M1 `0.2016`, M1c
`0.1693`, M2 `0.1342`, and M3
`0.1415`.

## Protein and target shuffle diagnostics

| diagnostic | median prediction correlation | median effect size (1-correlation) |
| --- | ---: | ---: |
| protein shuffle | `1.000000` | `0.000000` |
| target shuffle | `0.996871` | `0.003129` |

The protein and target shuffle files are inference-only and use no retraining
or model selection. The protein conditionality gate is **NO-GO** and the target
conditionality gate is **NO-GO** under the `0.995` threshold.

For direct historical comparison, v0.5 protein-shuffle median correlation was
`1.000000`
and target-shuffle median correlation was
`0.999871`.
The v0.6 values are in `v05_vs_v06_shuffle_comparison.csv`; both stages show
the same conditioning-collapse pattern.

## Seed stability and training health

Seed-level rows: 105 model/protein/seed results. Per-protein seed means and
standard deviations are in `seed_stability_per_protein.csv` and summary values
are in `seed_stability.csv`. All 60 model runs have positive prediction
variance and zero NaN/Inf values; detailed health is in `training_health.csv`.

## Gate table

```text
                            name  passed              observed                              threshold       stage
                  rc_class_count    True                  8192                                   8192  pre_replay
        v06_rc_class_injectivity    True          collisions=0            zero cross-class collisions  pre_replay
  synthetic_all_rankings_correct    True                   4/4                                    4/4  pre_replay
synthetic_protein_conditionality    True    min corr=-1.000000           positive correlation < 0.995  pre_replay
 synthetic_target_conditionality    True    min corr=-1.000000           positive correlation < 0.995  pre_replay
     held_out_spearman_available    True               105/105         105 complete protein-seed rows post_replay
     protein_shuffle_no_collapse   False  median corr=1.000000 median prediction correlation <= 0.995 post_replay
      target_shuffle_no_collapse   False  median corr=0.996871 median prediction correlation <= 0.995 post_replay
    prediction_variance_positive    True min variance=0.190398                      > 0 for every run post_replay
                      no_nan_inf    True         max NaN/Inf=0                                      0 post_replay
```

## Scientific interpretation

1. **Is the current failure mainly a representation bug?** Partly. The RC
   averaging bug is real and severe: 16,384 oriented 7-mers collapse to 2,000
   vectors, with 1,872 cross-class collision groups. Repairing it changes the
   ranking metrics, so it was not a harmless implementation detail. However,
   the post-repair shuffle collapse remains, so representation bug alone is not
   the complete explanation.
2. **How much did RC repair solve?** It restores exact RC-class injectivity and
   raises developmental replay median Spearman, but it does not restore robust
   protein conditioning.
3. **Is target representation important?** The old target representation is
   demonstrably information-losing; v0.6 preserves order/register/multiplicity.
   In this replay, target shuffle remains near-invariant, so importance for
   learned conditionality is not established.
4. **Does protein-conditioning collapse remain?** Yes. Protein shuffle median
   correlation is `1.000000`.
5. **Should sequence-based conditional modeling continue?** **STOP as a final
   model-development line.** The synthetic capability is present, but real
   protein/target effects are not used robustly and the sequence-only baseline
   remains competitive or better. A further small, pre-registered audit could
   test a new hypothesis, but no capacity/pair/attention scaling is justified.
6. **Next direction:** do not enter natural-to-designed transfer from this
   exposed cohort. The evidence supports a benchmark/generalization-gap paper
   or a carefully pre-registered structure-guided hypothesis after obtaining a
   new independent designed-DBP cohort. Natural-transfer H3 is already
   `NOT_SUPPORTED`; no confirmatory quantitative dataset was searched here.
