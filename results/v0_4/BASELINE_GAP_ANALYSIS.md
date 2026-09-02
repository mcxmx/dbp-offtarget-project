# v0.4 Baseline Gap Analysis

Analysis date: 2026-09-02

## Summary

The fixed benchmark is v0.3.1 GSE237017 designed-DBP uPBM: 7 proteins and 57,344 protein-RC-class experimental units. Scores are processed experimental PBM E-score consensus values and are evaluated as per-protein DNA ranking targets.

Best sequence-only baseline: `sequence_kmer3` with macro median Spearman 0.232. The empirical replicate agreement reference is median E-score replicate Spearman 0.591. The gap from the best sequence-only median to this reference is 0.359.

DeepPBS was not fairly evaluable in the current Windows environment because the official preprocessing stack requires additional structure-processing dependencies. SimpleProteinConditionalBaseline is intentionally untrained because no assay-matched natural PBM/uPBM training set has been added yet.

NA-MPNN ran only as a structural diagnostic for DBP35 and DBP48. Its macro median Spearman over these two evaluable proteins was -0.041. DBP48/8TAC appears in NA-MPNN split files and is not a zero-shot result.

## What Current Baselines Can Explain

- Sequence-only 3-mer similarity is the strongest Tier 0 baseline, but remains well below empirical replicate agreement.
- DBP35 NA-MPNN diagnostic Spearman is positive, but DBP48 is negative despite having an experimental structure and a training/validation overlap warning.
- Coverage is the main structure-aware bottleneck: five of seven designed DBPs lack a public structure/model found in the checked official sources.

## Disagreement Cases

v0.3.1 defined 1515 sequence-vs-experiment disagreement candidates. In v0.4, NA-MPNN predictions exist for 398 of these candidates and rank 50 of them in the top 10% of NA-MPNN scores. This is a diagnostic resolution rate of 0.126 among NA-MPNN-evaluable disagreement candidates, not across the whole benchmark.

The v0.4 failure table contains 1831 sequence-vs-experiment ranking cases: 1123 experimental-high/sequence-low cases and 708 sequence-high/experimental-low cases.

## Hardest Observed Regimes

For the best sequence-only baseline (`sequence_kmer3`), the lowest per-protein Spearman is -0.139 on DBP1. The lowest evaluated motif-distance stratum is DBP1 distance 3+ with Spearman -0.128 if using the same baseline.

## Per-Protein Metrics

| protein_id | baseline | spearman | evaluation_status |
| --- | --- | --- | --- |
| DBP1 | sequence_hamming | -0.194 | evaluated |
| DBP1 | sequence_edit | -0.194 | evaluated |
| DBP1 | sequence_kmer3 | -0.139 | evaluated |
| DBP1 | sequence_kmer4 | -0.072 | evaluated |
| DBP1 | NA-MPNN_structural_ppm | NA | not_evaluable_missing_public_structure |
| DBP1 | DeepPBS | NA | not_evaluable_current_environment |
| DBP1 | SimpleProteinConditionalBaseline | NA | not_trained_no_assay_matched_training_data |
| DBP3 | sequence_hamming | -0.024 | evaluated |
| DBP3 | sequence_edit | -0.027 | evaluated |
| DBP3 | sequence_kmer3 | -0.035 | evaluated |
| DBP3 | sequence_kmer4 | -0.018 | evaluated |
| DBP3 | NA-MPNN_structural_ppm | NA | not_evaluable_missing_public_structure |
| DBP3 | DeepPBS | NA | not_evaluable_current_environment |
| DBP3 | SimpleProteinConditionalBaseline | NA | not_trained_no_assay_matched_training_data |
| DBP35 | sequence_hamming | 0.011 | evaluated |
| DBP35 | sequence_edit | 0.019 | evaluated |
| DBP35 | sequence_kmer3 | 0.277 | evaluated |
| DBP35 | sequence_kmer4 | 0.154 | evaluated |
| DBP35 | NA-MPNN_structural_ppm | 0.250 | evaluated |
| DBP35 | DeepPBS | NA | not_evaluable_current_environment |
| DBP35 | SimpleProteinConditionalBaseline | NA | not_trained_no_assay_matched_training_data |
| DBP48 | sequence_hamming | 0.284 | evaluated |
| DBP48 | sequence_edit | 0.310 | evaluated |
| DBP48 | sequence_kmer3 | 0.345 | evaluated |
| DBP48 | sequence_kmer4 | 0.217 | evaluated |
| DBP48 | NA-MPNN_structural_ppm | -0.331 | evaluated |
| DBP48 | DeepPBS | NA | not_evaluable_current_environment |
| DBP48 | SimpleProteinConditionalBaseline | NA | not_trained_no_assay_matched_training_data |
| DBP5 | sequence_hamming | 0.084 | evaluated |
| DBP5 | sequence_edit | 0.089 | evaluated |
| DBP5 | sequence_kmer3 | 0.224 | evaluated |
| DBP5 | sequence_kmer4 | 0.100 | evaluated |
| DBP5 | NA-MPNN_structural_ppm | NA | not_evaluable_missing_public_structure |
| DBP5 | DeepPBS | NA | not_evaluable_current_environment |
| DBP5 | SimpleProteinConditionalBaseline | NA | not_trained_no_assay_matched_training_data |
| DBP6 | sequence_hamming | 0.056 | evaluated |
| DBP6 | sequence_edit | 0.066 | evaluated |
| DBP6 | sequence_kmer3 | 0.362 | evaluated |
| DBP6 | sequence_kmer4 | 0.288 | evaluated |
| DBP6 | NA-MPNN_structural_ppm | NA | not_evaluable_missing_public_structure |
| DBP6 | DeepPBS | NA | not_evaluable_current_environment |
| DBP6 | SimpleProteinConditionalBaseline | NA | not_trained_no_assay_matched_training_data |
| DBP9 | sequence_hamming | 0.067 | evaluated |
| DBP9 | sequence_edit | 0.061 | evaluated |
| DBP9 | sequence_kmer3 | 0.232 | evaluated |
| DBP9 | sequence_kmer4 | 0.232 | evaluated |
| DBP9 | NA-MPNN_structural_ppm | NA | not_evaluable_missing_public_structure |
| DBP9 | DeepPBS | NA | not_evaluable_current_environment |
| DBP9 | SimpleProteinConditionalBaseline | NA | not_trained_no_assay_matched_training_data |

## Macro Gap Table

| baseline | macro_median_spearman | macro_mean_spearman | n_proteins_evaluated | replicate_reference_spearman | gap_to_reference | disagreement_candidates_total |
| --- | --- | --- | --- | --- | --- | --- |
| DeepPBS | NA | NA | 0 | 0.591 | NA | 1515 |
| NA-MPNN_structural_ppm | -0.041 | -0.041 | 2 | 0.591 | 0.632 | 1515 |
| SimpleProteinConditionalBaseline | NA | NA | 0 | 0.591 | NA | 1515 |
| sequence_edit | 0.061 | 0.046 | 7 | 0.591 | 0.530 | 1515 |
| sequence_hamming | 0.056 | 0.041 | 7 | 0.591 | 0.536 | 1515 |
| sequence_kmer3 | 0.232 | 0.181 | 7 | 0.591 | 0.359 | 1515 |
| sequence_kmer4 | 0.154 | 0.129 | 7 | 0.591 | 0.438 | 1515 |

## Disagreement Resolution Table

| protein_id | baseline | n_v0_3_1_disagreement_candidates | n_evaluable_candidates | n_resolved | resolution_rate | resolution_definition | evaluation_status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DBP1 | NA-MPNN_structural_ppm | 287 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP1 | DeepPBS | 287 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP1 | SimpleProteinConditionalBaseline | 287 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP3 | NA-MPNN_structural_ppm | 381 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP3 | DeepPBS | 381 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP3 | SimpleProteinConditionalBaseline | 381 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP35 | NA-MPNN_structural_ppm | 174 | 174 | 50 | 0.287 | candidate predicted in top 10% within protein by protein-conditioned baseline score | evaluated |
| DBP35 | DeepPBS | 174 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP35 | SimpleProteinConditionalBaseline | 174 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP48 | NA-MPNN_structural_ppm | 224 | 224 | 0 | 0.000 | candidate predicted in top 10% within protein by protein-conditioned baseline score | evaluated |
| DBP48 | DeepPBS | 224 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP48 | SimpleProteinConditionalBaseline | 224 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP5 | NA-MPNN_structural_ppm | 178 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP5 | DeepPBS | 178 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP5 | SimpleProteinConditionalBaseline | 178 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP6 | NA-MPNN_structural_ppm | 150 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP6 | DeepPBS | 150 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP6 | SimpleProteinConditionalBaseline | 150 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP9 | NA-MPNN_structural_ppm | 121 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP9 | DeepPBS | 121 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |
| DBP9 | SimpleProteinConditionalBaseline | 121 | 0 | 0 | NA | candidate predicted in top 10% within protein by protein-conditioned baseline score | not_evaluable |

## Interpretation Limits

These results do not show that a new model is better than DeepPBS or NA-MPNN. They show that the current benchmark exposes a large sequence-only gap and that existing structure-aware methods are not yet comprehensively evaluable on the seven designed DBPs with public structures available in this repository.
