# Model Split Plan v2

Audit date: 2026-09-01

This v0.3.1 split plan supersedes the initial v0.3 draft for designed-DBP uPBM benchmark use.

## Split 1: Standard Protein-Level Split

Use only as a debugging split. Hold out proteins, but report that this does not guarantee independence if protein families, target groups, or assay types are shared.

## Split 2: Protein-Sequence-Cluster-Out

Cluster proteins by sequence identity and hold out entire clusters. For the seven GSE237017 designed DBPs, v0.3.1 uses a 0.60 identity threshold and identifies four clusters:

- DBP1/DBP3
- DBP5/DBP35
- DBP48
- DBP6/DBP9

DBP6 train and DBP9 test should not be treated as the strongest independent generalization test because they share high sequence identity and target/motif context.

## Split 3: Target-Group-Out

Hold out all proteins sharing an original target group, assay target group, or motif group. This guards against learning target-specific shortcuts rather than protein-conditioned recognition.

## Split 4: Designed DBP Zero-Shot External Test

Train on natural DBP experimental specificity data and evaluate DBP1, DBP3, DBP5, DBP6, DBP9, DBP35, and DBP48 as designed DBP external/OOD data.

This remains conditional on assay matching. If natural data use HT-SELEX and designed data use uPBM, the test measures protein shift plus assay shift.

## Split 5: Optional Leave-One-Designed-Protein-Out

Train/calibrate on six designed proteins and test one held-out designed protein. This is useful after a zero-shot benchmark is established, but it should not replace Split 4.

## DNA Split Rule

All DNA-level splits must group by canonical reverse-complement equivalence class.

Forbidden leakage:

- sequence in train, reverse complement in test
- canonical 7-mer class in both train and test
- motif-containing RC class split across train/test

## Primary Metrics

- Per-protein Spearman correlation against processed uPBM E-score.
- Per-protein top-k enrichment.
- Performance relative to empirical replicate agreement.
- OOD drop relative to assay-matched natural PBM/uPBM controls.
