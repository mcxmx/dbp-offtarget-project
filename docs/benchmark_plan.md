# Benchmark Plan

1. Positive: experimentally observed protein-DNA complexes from PDB, used as structural cognate pairs.
2. Negative: GC-matched random DNA and fully random DNA of identical length.
3. Hard negative: single mutants, double mutants, and genome near-matches.
4. Experimental ground truth sources: PDB, CIS-BP, PBM, HT-SELEX, and future protein-conditioned affinity datasets.
5. Leakage control: split by protein family, homology cluster, and DNA similarity, not by random row split.
6. Protein-family split: cluster proteins by sequence identity and keep families disjoint across train/test.
7. DNA-similarity split: cluster target and negative sequences by edit distance / k-mer similarity.
8. Sequence-only baseline: Hamming, edit distance, sequence identity, GC difference, k-mer similarity, reverse-complement-aware similarity.
9. Structure confidence comparison: add AF-derived or structure-derived confidence once available.
10. Genome-wide ranking: rank genomic loci by near-match distance, strand, and sequence proxy score.
11. Uncertainty validation: use model ensemble disagreement, seed disagreement, calibration curves, and coverage-risk tradeoffs once labels exist.

