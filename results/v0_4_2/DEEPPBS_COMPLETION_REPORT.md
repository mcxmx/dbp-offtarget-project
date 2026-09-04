# DeepPBS Baseline Completion Report

Audit date: 2026-09-04

## Scope

This report closes the previously incomplete DeepPBS baseline. The official
upstream checkout was run remotely on Ubuntu using the documented preprocessing
and prediction scripts. The project did not modify the upstream model or
weights. The PBM comparison is a per-protein ranking comparison between
processed experimental uPBM E-scores and a fixed PWM-derived sequence score.
The completion claim is technical and diagnostic: only structures with a
legitimate protein-DNA input are included, so coverage is reported explicitly.

## Structure coverage

| DBP | PDB | Evaluable | Reason | Overlap status |
|---|---|---:|---|---|
| DBP1 | NA | no | no reliable public protein-DNA complex structure/model | not evaluated; homolog audit unresolved |
| DBP3 | NA | no | no reliable public protein-DNA complex structure/model | not evaluated; homolog audit unresolved |
| DBP5 | NA | no | no reliable public protein-DNA complex structure/model | not evaluated; homolog audit unresolved |
| DBP6 | NA | no | no reliable public protein-DNA complex structure/model | not evaluated; homolog audit unresolved |
| DBP9 | NA | no | no reliable public protein-DNA complex structure/model | not evaluated; homolog audit unresolved |
| DBP35 | NA | yes | official run completed | no_exact_overlap_found_in_checked_manifests; homolog_unknown |
| DBP48 | 8TAC | yes | official run completed | no_exact_overlap_found_in_checked_manifests; homolog_unknown |

DeepPBS coverage is **2/7**. DBP35 uses the available theoretical
design complex. DBP48 uses the experimental 8TAC complex after a project-side
DSSR-defined helix-only input preparation; the original 8TAC file is preserved.

## Output semantics and scoring

The upstream source confirms `A/C/G/T` column order in
`deeppbs/dna_encodings.py::seqToOneHot` and `oneHotToSeq`. In
`run/predict.py`, `P` is the post-softmax ensemble output after averaging the
two strand halves with a reversed second half; `Seq` is the hard input
sequence. The project scorer evaluates every contiguous 7-mer PWM window with
`sum(log(P + 1e-9))`, takes the maximum window, then takes the maximum over a
candidate and its reverse complement. This is a fixed DeepPBS-derived ranking
proxy, not affinity, Kd, probability, or a calibrated binding score.

## PBM results

| DBP | Spearman | Status | RC classes |
|---|---:|---|---:|
| DBP1 | NA | not_evaluable_missing_prediction | NA |
| DBP3 | NA | not_evaluable_missing_prediction | NA |
| DBP5 | NA | not_evaluable_missing_prediction | NA |
| DBP6 | NA | not_evaluable_missing_prediction | NA |
| DBP9 | NA | not_evaluable_missing_prediction | NA |
| DBP35 | 0.040176 | evaluated | 8192 |
| DBP48 | 0.278366 | evaluated | 8192 |

DeepPBS evaluable macro median Spearman is **0.159271** across **2/7**
proteins. For context, the prior sequence-only k-mer3 median was
`0.232082`, SimpleProteinConditional was `0.361611`,
FrozenPLM was `0.153020`, and the empirical replicate agreement
reference was `0.591439`. These numbers have different coverage and
model assumptions; the two-protein DeepPBS median is diagnostic, not a
seven-protein generalization estimate.

## Hard-case integration

Using the frozen v0.3.1 disagreement definition, DeepPBS has
**398 eligible candidate rows**, resolves
**39**, leaves **359 unresolved**,
and has a resolution rate of **0.097990** within
its eligible denominator. The old 1,515 total is the all-protein candidate
set; it must not be used as DeepPBS's denominator because five proteins have
no legal structure input. The completed four-method common-hard subset contains
**25** rows and is written separately.

The frozen all-protein disagreement set contains
**1515** candidates. DeepPBS is not
assigned a denominator of 1,515 because DBP1, DBP3, DBP5, DBP6, and DBP9 have
no legal structure input. Candidate-level failure categories and pairwise
complementarity counts are in
`tables/baseline_failure_cases_deeppbs_completed_v0_4_2.parquet` and
`tables/baseline_complementarity_deeppbs_completed_v0_4_2.csv`.

## Baseline arena

The completed comparison table is
`tables/final_strong_baseline_summary_deeppbs_completed_v0_4_2.csv`.
DeepPBS is the only structure-aware result here and covers 2/7 proteins.
Its median Spearman is therefore a diagnostic statistic, not a fair
seven-protein ranking claim.

## Official example and run provenance

The official `5x6g` example passed in the Ubuntu VM. DBP35 completed with
Helix score 1.0 and contact count 266. DBP48 completed after DSSR-defined
helix-only preparation with Helix score 1.0 and contact count 224. Exact
commands, stdout/stderr, environment, input/output paths, and SHA256 values
are retained under `external_runs/` and in
`tables/deeppbs_run_manifest_completed_v0_4_2.csv`.

## Limitations

- Only DBP35 and DBP48 have a legal public/project structure input in this
  repository; five proteins remain not evaluable rather than being fabricated.
- DBP35 is a theoretical Rosetta model, not an experimental complex.
- DBP48 required a transparent helix-only preprocessing repair because the
  deposited complex has non-helical DNA overhang residues that break the
  upstream shape extraction path.
- The checked DeepPBS manifests did not contain exact designed-DBP/PDB hits,
  but full homolog-level training-set audit remains unresolved because the
  upstream repository does not distribute all training sequences.
- A two-protein result cannot establish a systematic strong-baseline ranking.
