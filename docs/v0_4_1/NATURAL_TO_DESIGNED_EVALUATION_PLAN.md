# v0.4.1 Natural-to-Designed Evaluation Plan

Date: 2026-09-02

## Train/Test Design

- Train: UniPROBE natural PBM train proteins from `metadata/v0_4_1/natural_pbm_splits.csv`.
- Validation: UniPROBE natural PBM validation proteins; no designed DBP rows used for hyperparameter choice.
- Test A: held-out natural PBM proteins.
- Test B: GSE237017 designed DBPs, kept external.

## Confounders

- Natural UniPROBE uses processed contiguous 8-mer E-scores.
- Designed GSE237017 uses processed uPBM 7-mer E-scores.
- Protein sequences in v0.4.1 are full-length UniProt references, not confirmed assay constructs.
- A natural-to-designed performance drop may reflect assay/platform/k-mer processing shift as well as protein-distribution shift.

## Split Rules

- Protein split is cluster-aware at the 40% proxy identity level.
- DNA rows are grouped by reverse-complement canonical class; no oriented RC pair is split into separate evaluation units.
- Designed DBP sequences are excluded from natural training.
