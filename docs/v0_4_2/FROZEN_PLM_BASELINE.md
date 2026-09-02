# v0.4.2 FrozenPLMProteinConditionalBaseline

Date: 2026-09-02

This baseline uses frozen ESM-2 `esm2_t12_35M_UR50D` mean-pooled protein embeddings concatenated with simple DNA features and scored by ridge regression. It is a baseline only, not the proposed model.

## Training and validation

- Training data: natural PBM train proteins only.
- Validation: natural PBM validation proteins only.
- External test: GSE237017 designed DBPs only after alpha selection.
- Selected alpha: 1000.0
- Validation macro median Spearman: 0.287

The protein LM is not fine-tuned. Designed DBP rows do not influence checkpoint, alpha, feature encoding, or any hyperparameter selection.
