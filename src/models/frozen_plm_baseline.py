from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.models.base_binding_model import BaseBindingModel
from src.models.simple_protein_conditional_baseline import dna_kmer_features


@dataclass
class FrozenPLMProteinConditionalBaseline(BaseBindingModel):
    """Frozen protein language model baseline, not a proposed architecture.

    The protein representation is computed upstream by a frozen ESM-2 model and
    is never fine-tuned on either natural PBM or designed uPBM data. The scoring
    head is intentionally small and deterministic so it can serve as a strong
    baseline closure experiment rather than a new method.
    """

    is_protein_conditioned: bool = True
    protein_lm_frozen: bool = True
    not_proposed_method: bool = True
    checkpoint_name: str = "esm2_t12_35M_UR50D"
    score_label: str = "frozen_plm_score"
    weights: np.ndarray | None = None
    intercept: float = 0.0
    feature_mean: np.ndarray | None = None
    feature_scale: np.ndarray | None = None
    protein_projection: np.ndarray | None = None

    def featurize(self, protein_embedding: np.ndarray, dna_sequence: str) -> np.ndarray:
        if self.protein_projection is None:
            raise RuntimeError("protein_projection is required for FrozenPLMProteinConditionalBaseline")
        projected = protein_embedding @ self.protein_projection
        dna = dna_kmer_features(dna_sequence)
        interaction = np.outer(projected, dna).ravel()
        return np.concatenate([projected, dna, interaction]).astype(np.float32)

    def score(self, protein_embedding: np.ndarray, dna_sequence: str) -> float:
        if self.weights is None or self.feature_mean is None or self.feature_scale is None:
            raise RuntimeError("FrozenPLMProteinConditionalBaseline head is not fitted")
        features = self.featurize(protein_embedding, dna_sequence)
        scaled = (features - self.feature_mean) / self.feature_scale
        return float(scaled @ self.weights + self.intercept)
