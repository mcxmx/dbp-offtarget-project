# Benchmark Quality Report v0.2

Generated from repository tables on 2026-08-29.

## Summary Counts

- Raw historical PDB pairs: 16
- Curated v0.2 core benchmark pairs: 8
- Sequence-specific pairs: 7
- Non-specific pairs: 2
- Guide-dependent pairs: 2
- Lesion-specific pairs: 2
- Designed sequence-specific pairs: 2
- Uncertain pairs: 0
- PDB pairs with quantitative specificity ground truth: 0
- PDB pairs without quantitative specificity ground truth: 16
- Experimental specificity pilot proteins: 5
- Experimental specificity pilot rows: 1209

## Sequence Specificity Classes

| sequence_specificity_class | count |
| --- | --- |
| sequence_specific | 7 |
| designed_sequence_specific | 2 |
| lesion_specific | 2 |
| guide_dependent | 2 |
| non_specific | 2 |
| structure_specific | 1 |

## Recommended Use

| recommended_use | count |
| --- | --- |
| core_benchmark | 8 |
| exclude_from_specificity_benchmark | 4 |
| negative_control | 2 |
| auxiliary_case | 1 |
| method_demo_only | 1 |

## Benchmark v0.2 Candidate Counts

| candidate_type | count |
| --- | --- |
| double_mut | 9360 |
| gc_matched_random | 8000 |
| random_dna | 8000 |
| single_mut | 393 |
| target | 8 |
| total | 25761 |

## v0.1 to v0.2 Corrections

- Split PDB structural evidence from quantitative specificity ground truth.
- Set all PDB-only records to `has_quantitative_specificity_ground_truth=False`.
- Replaced the longest-chain benchmark assumption with curated chain annotation
  plus chain-contact evidence where available.
- Moved guide-dependent, lesion-specific, non-specific, and transposase/substrate
  cases out of the v0.2 core specificity benchmark.
- Replotted v0.2 figures with sequence-only proxy terminology.
- Added an explicit positional-bias analysis for k-mer and combined proxy
  metrics.

## Current Interpretation

The v0.2 structural/mutation benchmark is suitable for reproducible pipeline
testing and sequence-only sanity checks. It is not yet a calibrated
protein-conditioned off-target predictor. The experimental specificity pilot is
a separate Layer C resource based on JASPAR PFM-derived PWM scores, not raw PBM
or HT-SELEX enrichment.
