# Benchmark Plan

1. Positive: PDB-derived structural cognate pairs for Layer A, plus experimentally measured specificity records for Layer C when quantitative assay data are available.
2. Negative: GC-matched random DNA and fully random DNA of identical length.
3. Hard negative: single mutants, double mutants, and genome near-matches.
4. Experimental ground truth sources: PBM, HT-SELEX, CIS-BP/JASPAR-like curated motif resources, Kd/competition assays, and future affinity datasets. PDB alone is structural evidence, not quantitative specificity ground truth.
5. Leakage control: split by protein family, homology cluster, and DNA similarity, not by random row split.
6. Protein-family split: cluster proteins by sequence identity and keep families disjoint across train/test.
7. DNA-similarity split: cluster target and negative sequences by edit distance / k-mer similarity.
8. Sequence-only baseline: Hamming, edit distance, sequence identity, GC difference, k-mer similarity, reverse-complement-aware similarity.
9. Structure confidence comparison: add AF-derived or structure-derived confidence once available.
10. Genome-wide ranking: retrieve genomic candidate loci by near-match distance and strand first; call it off-target prediction only after protein-conditioned scoring or experimental validation is added.
11. Uncertainty validation: use model ensemble disagreement, seed disagreement, calibration curves, and coverage-risk tradeoffs once labels exist.
