# v0.4.1 Natural PBM QC Rules

Audit date: 2026-09-02

## Unit Definitions

- Natural PBM source rows are UniPROBE processed contiguous 8-mer E-score rows.
- The independent DNA unit used for splitting/evaluation is the reverse-complement canonical 8-mer class.
- A complete contiguous 8-mer profile is expected to contain 32,896 RC classes.

## Quality Levels

- `high_confidence`: complete RC-class coverage, DNA/RC QC pass, clear non-fusion protein label, and high-confidence UniProt reference sequence.
- `usable`: complete RC-class coverage, DNA/RC QC pass, clear non-fusion protein label, and medium-confidence UniProt reference sequence.
- `low_confidence`: complete/incomplete PBM profile that lacks a conservative protein sequence, or incomplete 8-mer coverage.
- `exclude`: invalid DNA/RC QC or complex/fusion/unclear labels.

Only `high_confidence` and `usable` rows enter `natural_pbm_benchmark_v0_4_1.parquet`.
Full-length UniProt sequences are retained with `sequence_match_to_assay=false` unless construct-level sequence evidence is available.

Protein clusters in v0.4.1 are generated from a fast amino-acid 3-mer Jaccard proxy to keep obvious duplicate and near-duplicate reference sequences in the same split. This is a split-hygiene proxy, not a replacement for a future MMseqs2/CD-HIT homology audit.
