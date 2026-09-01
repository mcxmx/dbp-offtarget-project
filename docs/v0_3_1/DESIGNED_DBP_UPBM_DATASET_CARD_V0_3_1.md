# Designed DBP uPBM Dataset Card v0.3.1

Dataset: GSE237017 designed-DBP processed experimental uPBM benchmark.

Audit date: 2026-09-01

## Source

- Paper: "Computational design of sequence-specific DNA-binding proteins"
- DOI: 10.1038/s41594-025-01669-4
- GEO: GSE237017
- Source-data files: Source Data Fig. 4 and Source Data Extended Data Fig. 8
- Supplementary tables: official Nature supplementary workbook

## Units

- Proteins: 7
- GEO usable samples: 12
- Oriented 7-mer rows per protein: 16384
- Reverse-complement equivalence classes per protein: 8192
- Total oriented protein-7mer rows: 114688
- Total protein-RC-class experimental units: 57344

The 114688 oriented rows are not 114688 independent DNA sequence units. They collapse to 57344 protein-RC-class units because uPBM measures dsDNA and the processed files pair each 7-mer with its reverse complement.

## Primary Score

Primary score in v0.3.1: `experimental_escore_consensus`.

Definition: replicate consensus mean of GEO processed uPBM E-score for a protein and 7-mer/RC class.

Not justified:

- Kd
- binding free energy
- binding probability
- in vivo genomic binding
- absolute cross-protein affinity

## Target Definitions

v0.3.1 separates:

- original design target
- experimental assay target
- PBM evaluation motif

DBP48 is the important correction:

- Original design target: target I (`CGCCCAAAGCCGCG`)
- Experimental assay/PBM evaluation reference: target C (`CGACACCTGACGCG`)
- PBM motif used for reproduction: `CTGACG`

See `metadata/v0_3_1/designed_dbp_target_definitions.csv`.

## Paper Percentile Reproduction

Extended Data Fig. 8 motif percentiles were reproduced using 8192 source-data RC-class rows and a <=2 percentile point tolerance.

Result: PASS for 7/7 proteins.

Maximum absolute difference: 1.2535 percentile points.

See `results/v0_3_1/tables/paper_percentile_reproduction.csv`.

## RC-Aware Baseline

v0.3.1 computes sequence-only metrics against the paper PBM motif while considering candidate and motif reverse-complement equivalence. Metrics are computed per protein on 8192 RC classes.

This remains a sequence-only proxy baseline and is not protein-conditioned.

## Independence

Protein sequence clusters at identity threshold 0.60:

- protein_cluster_1: DBP1, DBP3
- protein_cluster_2: DBP5, DBP35
- protein_cluster_3: DBP48
- protein_cluster_4: DBP6, DBP9

Shared target/motif groups must be respected in future splits.

## Recommended Use

- External designed-DBP benchmark for v0.4 protein-conditioned baselines.
- Per-protein ranking against processed uPBM E-score.
- Assay reproducibility-aware model comparison using replicate agreement.

## Prohibited Interpretation

- Do not claim full-target affinity from overlapping 7-mer or motif summaries.
- Do not pool all proteins and report a single absolute-score regression metric.
- Do not split reverse-complement equivalent DNA sequences across train/test.
- Do not interpret natural-to-designed performance drop without assay-shift controls.
