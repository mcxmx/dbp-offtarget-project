from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from statistics import mean, pstdev
from typing import Any, Iterable

from src.utils import canonical_kmer_set, clamp, gc_content


@dataclass
class PredictionBundle:
    mean_score: float | None = None
    std_score: float | None = None
    coefficient_of_variation: float | None = None
    model_disagreement: float | None = None
    seed_disagreement: float | None = None
    entropy: float | None = None
    scores: list[float] = field(default_factory=list)
    details: dict[str, Any] = field(default_factory=dict)


class BaseBindingModel(ABC):
    is_protein_conditioned: bool | None = None
    score_label: str = "unspecified_model_score"

    @abstractmethod
    def score(self, protein_sequence: str, dna_sequence: str) -> float:
        raise NotImplementedError

    def predict_with_uncertainty(
        self,
        protein_sequence: str,
        dna_sequence: str,
        models: Iterable["BaseBindingModel"] | None = None,
        seeds: Iterable[int] | None = None,
    ) -> PredictionBundle:
        scores: list[float] = []
        if models is not None:
            for model in models:
                scores.append(model.score(protein_sequence, dna_sequence))
        else:
            scores.append(self.score(protein_sequence, dna_sequence))
        if seeds is not None:
            for _ in seeds:
                scores.append(self.score(protein_sequence, dna_sequence))
        if not scores:
            return PredictionBundle()
        mu = mean(scores)
        sigma = pstdev(scores) if len(scores) > 1 else 0.0
        cv = sigma / mu if mu not in (0.0, None) else None
        return PredictionBundle(
            mean_score=mu,
            std_score=sigma,
            coefficient_of_variation=cv,
            model_disagreement=sigma,
            seed_disagreement=sigma,
            entropy=None,
            scores=scores,
            details={"n_scores": len(scores), "is_protein_conditioned": self.is_protein_conditioned},
        )


class SequenceProxyBaseline(BaseBindingModel):
    """DNA-sequence-only sanity-check baseline; not a protein-conditioned model."""

    is_protein_conditioned = False
    score_label = "sequence_only_proxy_score"

    def __init__(self, k_values: tuple[int, ...] = (3, 4), weights: dict[str, float] | None = None):
        self.k_values = tuple(k_values)
        self.weights = weights or {
            "identity": 0.40,
            "edit_similarity": 0.15,
            "kmer3": 0.20,
            "kmer4": 0.15,
            "rc_kmer4": 0.10,
        }

    def score(self, protein_sequence: str, dna_sequence: str) -> float:
        dna_sequence = dna_sequence.upper().replace(" ", "")
        if not dna_sequence:
            return 0.0
        gc_term = 1.0 - abs(gc_content(dna_sequence) - 0.5)
        k4 = len(canonical_kmer_set(dna_sequence, 4, rc_aware=True)) / max(len(dna_sequence) - 3, 1)
        length_term = min(len(dna_sequence), 50) / 50.0
        score = 0.45 * gc_term + 0.35 * k4 + 0.20 * length_term
        return clamp(score, 0.0, 1.0)
