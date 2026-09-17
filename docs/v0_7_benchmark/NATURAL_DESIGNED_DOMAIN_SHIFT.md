# Natural–Designed Domain Shift

This is a data-layer comparison, not a new model or transfer run. Natural PBM uses UniPROBE contiguous 8-mers (57 proteins); designed data uses GSE237017 uPBM 7-mers (7 proteins). The local audit reports no verified exact assay constructs for the natural proteins, so protein sequence comparisons are reference-sequence diagnostics.

| metric | dataset | value_median | value_mean | n_rows |
| --- | --- | --- | --- | --- |
| aa_fraction_A | designed | 0.1277 | 0.1277 | 1 |
| aa_fraction_A | natural | 0.0680 | 0.0680 | 1 |
| aa_fraction_C | designed | 0.0026 | 0.0026 | 1 |
| aa_fraction_C | natural | 0.0155 | 0.0155 | 1 |
| aa_fraction_D | designed | 0.0174 | 0.0174 | 1 |
| aa_fraction_D | natural | 0.0524 | 0.0524 | 1 |
| aa_fraction_E | designed | 0.1252 | 0.1252 | 1 |
| aa_fraction_E | natural | 0.0573 | 0.0573 | 1 |
| aa_fraction_F | designed | 0.0091 | 0.0091 | 1 |
| aa_fraction_F | natural | 0.0334 | 0.0334 | 1 |
| aa_fraction_G | designed | 0.0709 | 0.0709 | 1 |
| aa_fraction_G | natural | 0.0535 | 0.0535 | 1 |
| aa_fraction_H | designed | 0.0025 | 0.0025 | 1 |
| aa_fraction_H | natural | 0.0294 | 0.0294 | 1 |
| aa_fraction_I | designed | 0.0661 | 0.0661 | 1 |
| aa_fraction_I | natural | 0.0442 | 0.0442 | 1 |
| aa_fraction_K | designed | 0.0740 | 0.0740 | 1 |
| aa_fraction_K | natural | 0.0598 | 0.0598 | 1 |
| aa_fraction_L | designed | 0.1054 | 0.1054 | 1 |
| aa_fraction_L | natural | 0.0823 | 0.0823 | 1 |
| aa_fraction_M | designed | 0.0022 | 0.0022 | 1 |
| aa_fraction_M | natural | 0.0227 | 0.0227 | 1 |
| aa_fraction_N | designed | 0.0194 | 0.0194 | 1 |
| aa_fraction_N | natural | 0.0666 | 0.0666 | 1 |
| aa_fraction_P | designed | 0.0140 | 0.0140 | 1 |
| aa_fraction_P | natural | 0.0676 | 0.0676 | 1 |
| aa_fraction_Q | designed | 0.0284 | 0.0284 | 1 |
| aa_fraction_Q | natural | 0.0527 | 0.0527 | 1 |
| aa_fraction_R | designed | 0.1378 | 0.1378 | 1 |
| aa_fraction_R | natural | 0.0491 | 0.0491 | 1 |
| aa_fraction_S | designed | 0.0481 | 0.0481 | 1 |
| aa_fraction_S | natural | 0.1080 | 0.1080 | 1 |
| aa_fraction_T | designed | 0.0567 | 0.0567 | 1 |
| aa_fraction_T | natural | 0.0631 | 0.0631 | 1 |
| aa_fraction_V | designed | 0.0703 | 0.0703 | 1 |
| aa_fraction_V | natural | 0.0420 | 0.0420 | 1 |
| aa_fraction_W | designed | 0.0076 | 0.0076 | 1 |
| aa_fraction_W | natural | 0.0058 | 0.0058 | 1 |
| aa_fraction_Y | designed | 0.0147 | 0.0147 | 1 |
| aa_fraction_Y | natural | 0.0267 | 0.0267 | 1 |
| dna_kmer_length | designed | 7.0000 | 7.0000 | 7 |
| dna_kmer_length | natural | 8.0000 | 8.0000 | 57 |
| nearest_natural_embedding_cosine_distance | designed | 0.2220 | 0.2410 | 7 |
| nearest_natural_embedding_euclidean_distance | designed | 3.9172 | 3.9702 | 7 |
| nearest_natural_sequence_similarity | designed | 0.1687 | 0.1676 | 7 |
| protein_sequence_length | designed | 58.0000 | 59.8571 | 7 |
| protein_sequence_length | natural | 505.0000 | 555.7193 | 57 |
| score_95_5_spread | designed | 0.3615 | 0.3773 | 7 |
| score_95_5_spread | natural | 0.4606 | 2.5856 | 57 |
| score_dynamic_range | designed | 0.8018 | 0.8194 | 7 |
| score_dynamic_range | natural | 0.9490 | 563.2783 | 57 |
| top1pct_position_entropy_mean | designed | 1.7559 | 1.7612 | 7 |
| top1pct_position_entropy_mean | natural | 1.7959 | 1.7759 | 57 |

## Interpretation

The two datasets differ in k-mer length, assay source, score processing and construct provenance. Frozen ESM-2 embeddings and rough sequence-similarity nearest neighbors are provided to test whether designed proteins occupy the natural reference distribution, but no PCA/UMAP distance is treated as a formal domain-shift test. Score dynamic range and entropy are kept assay-specific and are not pooled across assays. The correct conclusion is a plausible, measurable domain/construct/assay shift hypothesis, not proof that domain shift alone causes the generalization gap.
