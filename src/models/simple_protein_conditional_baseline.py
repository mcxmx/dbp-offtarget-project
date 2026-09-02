from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from src.models.base_binding_model import BaseBindingModel
from src.utils import DNA_ALPHABET, PROTEIN_ALPHABET, normalize_sequence


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
DINUC_ORDER = [a + b for a in DNA_ALPHABET for b in DNA_ALPHABET]
MAX_DNA_LEN = 8


def protein_composition_features(protein_sequence: str) -> np.ndarray:
    """Small fixed protein representation for non-final baseline prototypes."""
    seq = normalize_sequence(protein_sequence)
    length = max(len(seq), 1)
    counts = np.array([seq.count(aa) / length for aa in AA_ORDER], dtype=float)
    charged = (seq.count("D") + seq.count("E") + seq.count("K") + seq.count("R")) / length
    basic = (seq.count("K") + seq.count("R") + seq.count("H")) / length
    acidic = (seq.count("D") + seq.count("E")) / length
    valid_fraction = sum(residue in PROTEIN_ALPHABET for residue in seq) / length
    return np.concatenate([counts, np.array([length / 100.0, charged, basic, acidic, valid_fraction])])


def dna_kmer_features(dna_sequence: str, max_len: int = MAX_DNA_LEN) -> np.ndarray:
    seq = normalize_sequence(dna_sequence)
    if len(seq) < 1 or len(seq) > max_len or any(base not in DNA_ALPHABET for base in seq):
        raise ValueError(
            f"SimpleProteinConditionalBaseline expects an ACGT sequence with length 1-{max_len}, got {dna_sequence!r}"
        )
    one_hot = np.zeros((max_len, 4), dtype=float)
    base_to_idx = {base: i for i, base in enumerate(DNA_ALPHABET)}
    for i, base in enumerate(seq):
        one_hot[i, base_to_idx[base]] = 1.0
    dinuc = np.zeros(len(DINUC_ORDER), dtype=float)
    dinuc_index = {dinuc_value: i for i, dinuc_value in enumerate(DINUC_ORDER)}
    for i in range(len(seq) - 1):
        dinuc[dinuc_index[seq[i : i + 2]]] += 1.0
    dinuc /= max(len(seq) - 1, 1)
    base_comp = np.array([seq.count(base) / len(seq) for base in DNA_ALPHABET], dtype=float)
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    return np.concatenate([one_hot.ravel(), base_comp, dinuc, np.array([gc, len(seq) / max_len])])


def dna_7mer_features(dna_sequence: str) -> np.ndarray:
    seq = normalize_sequence(dna_sequence)
    if len(seq) != 7:
        raise ValueError(f"dna_7mer_features expects a 7-mer, got {dna_sequence!r}")
    return dna_kmer_features(seq)


def simple_pc_features(protein_sequence: str, dna_sequence: str) -> np.ndarray:
    protein = protein_composition_features(protein_sequence)
    dna = dna_kmer_features(dna_sequence)
    # A small interaction term keeps the baseline protein-conditioned without
    # becoming a proposed architecture.
    interaction = np.array([
        protein[20] * dna[-1],
        protein[22] * dna[-1],
        protein[23] * dna[-1],
    ])
    return np.concatenate([protein, dna, interaction])


@dataclass
class SimpleProteinConditionalBaseline(BaseBindingModel):
    """Minimal protein-conditioned baseline scaffold; not the proposed model.

    v0.4 intentionally does not train this model on GSE237017 as a main result,
    because no assay-compatible natural PBM/uPBM training set is in the project
    yet. The class exists to make the future Tier 1 interface explicit.
    """

    is_protein_conditioned: bool = True
    score_label: str = "simple_protein_conditional_baseline_score"
    weights: np.ndarray | None = None
    intercept: float = 0.0
    training_status: str = "untrained_no_assay_matched_natural_training_data"
    feature_names: list[str] = field(default_factory=list)

    def score(self, protein_sequence: str, dna_sequence: str) -> float:
        if self.weights is None:
            raise RuntimeError(
                "SimpleProteinConditionalBaseline is untrained. v0.4 records this "
                "as a feasibility scaffold, not an evaluated baseline."
            )
        features = simple_pc_features(protein_sequence, dna_sequence)
        if features.shape[0] != self.weights.shape[0]:
            raise ValueError("Feature vector and weight vector have different lengths")
        return float(features @ self.weights + self.intercept)
