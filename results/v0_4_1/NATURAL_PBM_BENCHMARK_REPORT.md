# v0.4.1 Natural PBM Benchmark Report

Audit date: 2026-09-02

## Summary

- Source: UniPROBE processed contiguous 8-mer E-score profiles.
- Final train/evaluation benchmark proteins: 57
- Protein families/coarse classes: 7
- Species represented: 7
- Protein-DNA units: 1,875,072 protein-RC-class rows.
- DNA unit: 8-mer reverse-complement equivalence class.
- Score: UniPROBE contiguous 8-mer E-score; used for per-protein ranking, not cross-protein absolute affinity.
- Protein sequence completeness in final benchmark: 100.0%
- Construct sequence completeness: 0.0%; reference full-length sequences are used and flagged as not assay-construct matched.

## QC Counts

quality_level
high_confidence    49
low_confidence     25
exclude            16
usable              8

## Cluster-Aware Split

       split  n_proteins  n_families  n_protein_dna_units
natural_test           9           5               296064
       train          39           7              1282944
  validation           9           4               296064

## Training Readiness

The benchmark is sufficient for a first simple protein-conditioned baseline because it has >50 proteins with complete 8-mer profiles and conservative reference sequences. It is not yet sufficient for claims about assay-free biological OOD because natural UniPROBE 8-mer profiles and designed GSE237017 uPBM 7-mer profiles differ in k-mer length, protocol, and score processing.
