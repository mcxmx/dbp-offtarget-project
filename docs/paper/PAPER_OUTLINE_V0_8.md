# Paper Outline V0.8

## Introduction
Define transfer of DNA-binding specificity to unseen de novo proteins and distinguish performance from conditionality.

## Result 1: Strict unseen-protein benchmark
- Claim: current sequence-based models show modest transfer and a reproducible benchmark gap.
- Supporting files: `results/v0_7_benchmark/benchmark_master.csv`, `results/v0_8_sota/sota_benchmark_master.csv`.
- Figure: Figure 1, cohort/split schematic.
- Limitation: seven exposed designed proteins.
- Reviewer objection: protein-level N is small and methods are heterogeneous.

## Result 2: Representation repair
- Claim: RC collision repair improves ranking but does not restore conditionality.
- Supporting files: `docs/v0_6_representation/`, `results/v0_6_representation/`.
- Figure: Figure 2, performance vs conditionality.
- Limitation: replay remains development-exposed.
- Reviewer objection: repair is necessary but not sufficient.

## Result 3: Failed rescue attempts
- Claim: dense supervision and local/residue attention did not rescue conditioning collapse.
- Supporting files: `results/v0_5_dense/`, `results/v0_5_local/`.
- Figure: Figure 3.
- Limitation: pilots are not universal architecture proof.
- Reviewer objection: other structure-aware methods were not yet fully covered.

## Result 4: Structure-method audit
- Claim: DeepPBS/NA-MPNN are legal partial diagnostics but current coverage is N=2; AF3 contact metrics are feasibility-only.
- Supporting files: `docs/v0_8_sota/`, `results/v0_8_sota/`.
- Figure: coverage-aware benchmark matrix.
- Limitation: cannot infer universal structure-method failure.
- Reviewer objection: request broader structure coverage.

## Result 5: Hard cases and domain shift
- Claim: 871/1515 common sequence-model failures and measurable natural/designed assay/construct differences motivate, but do not prove, a domain-shift/data-regime explanation.
- Supporting files: v0.7 hard taxonomy and shift CSV.
- Figure: Figure 4.
- Limitation: causal contributions are unidentified.
- Reviewer objection: hard cases are label/noise dependent.

## Result 6: Independent validation
- Status: blocked; do not include a Figure 5 until a preregistered cohort is available.

## Discussion / Limitations / Methods
Emphasize protein-level N=7, development exposure, no universal computational claim, and the highest-value next experiment: an independent assay-matched cohort with >7 designed proteins and dense specificity landscapes.
