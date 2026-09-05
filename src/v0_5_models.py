from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import torch
from torch import nn

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.utils import edit_distance, kmer_jaccard, normalize_sequence


DNA_ALPHABET = "ACGT"
DNA_INDEX = {base: index for index, base in enumerate(DNA_ALPHABET)}
KMER_LENGTH = 7


def validate_dna_sequence(sequence: str, *, min_length: int = 1) -> str:
    normalized = normalize_sequence(sequence)
    if len(normalized) < min_length or any(base not in DNA_INDEX for base in normalized):
        raise ValueError(f"Expected ACGT DNA of length >= {min_length}, got {sequence!r}")
    return normalized


def target_windows(target: str, k: int = KMER_LENGTH) -> tuple[str, ...]:
    """Return unique canonical RC windows from an orientation-agnostic target."""
    target = validate_dna_sequence(target, min_length=k)
    windows = {canonical_rc(target[index : index + k]) for index in range(len(target) - k + 1)}
    if not windows:
        raise ValueError(f"Target must contain at least one {k}-mer window")
    return tuple(sorted(windows))


def _one_hot(sequence: str) -> np.ndarray:
    sequence = validate_dna_sequence(sequence)
    result = np.zeros((len(sequence), len(DNA_ALPHABET)), dtype=np.float32)
    for index, base in enumerate(sequence):
        result[index, DNA_INDEX[base]] = 1.0
    return result.reshape(-1)


def rc_symmetric_one_hot(sequence: str) -> np.ndarray:
    """Encode both orientations and average them to make orientation irrelevant."""
    sequence = validate_dna_sequence(sequence)
    return ((_one_hot(sequence) + _one_hot(reverse_complement(sequence))) / 2.0).astype(np.float32)


def batch_rc_symmetric_one_hot(sequences: Sequence[str]) -> torch.Tensor:
    if not sequences:
        raise ValueError("At least one DNA sequence is required")
    lengths = {len(validate_dna_sequence(sequence)) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("A batch must contain equal-length DNA sequences")
    return torch.from_numpy(np.vstack([rc_symmetric_one_hot(sequence) for sequence in sequences]))


def target_window_features(target: str, k: int = KMER_LENGTH) -> torch.Tensor:
    return batch_rc_symmetric_one_hot(target_windows(target, k=k))


class SymmetricDNAEncoder(nn.Module):
    """Small shared DNA projection; RC symmetrization occurs before projection."""

    def __init__(self, input_dim: int = KMER_LENGTH * 4, hidden_dim: int = 32) -> None:
        super().__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.projection = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
        )

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        if features.ndim not in (2, 3) or features.shape[-1] != self.input_dim:
            raise ValueError(
                f"Expected DNA features [batch,{self.input_dim}] or "
                f"[batch,windows,{self.input_dim}], got {tuple(features.shape)}"
            )
        original_shape = features.shape
        projected = self.projection(features.reshape(-1, self.input_dim))
        return projected.reshape(*original_shape[:-1], self.hidden_dim)


class ProteinProjection(nn.Module):
    def __init__(self, protein_dim: int, hidden_dim: int) -> None:
        super().__init__()
        self.projection = nn.Sequential(
            nn.Linear(protein_dim, hidden_dim),
            nn.Tanh(),
        )

    def forward(self, protein_embedding: torch.Tensor) -> torch.Tensor:
        return self.projection(protein_embedding)


def _mlp(input_dim: int, hidden_dim: int, output_dim: int = 1) -> nn.Sequential:
    return nn.Sequential(
        nn.Linear(input_dim, hidden_dim),
        nn.ReLU(),
        nn.Linear(hidden_dim, output_dim),
    )


class CandidateDNAOnly(nn.Module):
    model_name = "M0_CandidateDNAOnly"
    inputs = "D"

    def __init__(self, hidden_dim: int = 32) -> None:
        super().__init__()
        self.dna_encoder = SymmetricDNAEncoder(hidden_dim=hidden_dim)
        self.head = _mlp(hidden_dim, hidden_dim)

    def forward(self, candidate_features: torch.Tensor) -> torch.Tensor:
        return self.head(self.dna_encoder(candidate_features)).squeeze(-1)


class ProteinCandidate(nn.Module):
    model_name = "M1_ProteinCandidate"
    inputs = "P,D"

    def __init__(self, protein_dim: int, hidden_dim: int = 32) -> None:
        super().__init__()
        self.protein_projection = ProteinProjection(protein_dim, hidden_dim)
        self.dna_encoder = SymmetricDNAEncoder(hidden_dim=hidden_dim)
        self.head = _mlp(hidden_dim * 3, hidden_dim)

    def forward(
        self,
        protein_embedding: torch.Tensor,
        candidate_features: torch.Tensor,
    ) -> torch.Tensor:
        protein = self.protein_projection(protein_embedding)
        candidate = self.dna_encoder(candidate_features)
        interaction = torch.cat([protein, candidate, protein * candidate], dim=-1)
        return self.head(interaction).squeeze(-1)


class TargetCandidateOnly(nn.Module):
    model_name = "M2_TargetCandidateOnly"
    inputs = "T,D"

    def __init__(self, hidden_dim: int = 32) -> None:
        super().__init__()
        self.dna_encoder = SymmetricDNAEncoder(hidden_dim=hidden_dim)
        self.comparative = _mlp(hidden_dim * 4, hidden_dim, hidden_dim)
        self.head = nn.Linear(hidden_dim, 1)

    def forward(
        self,
        target_window_features_batch: torch.Tensor,
        candidate_features: torch.Tensor,
    ) -> torch.Tensor:
        target = self.dna_encoder(target_window_features_batch).mean(dim=1)
        candidate = self.dna_encoder(candidate_features)
        comparative = torch.cat(
            [candidate, target, candidate - target, candidate * target],
            dim=-1,
        )
        return self.head(self.comparative(comparative)).squeeze(-1)


class ProteinTargetCandidate(nn.Module):
    """Minimal non-separable P,T,D model with protein-controlled FiLM."""

    model_name = "M3_ProteinTargetCandidate"
    inputs = "P,T,D"

    def __init__(self, protein_dim: int, hidden_dim: int = 32) -> None:
        super().__init__()
        self.protein_projection = ProteinProjection(protein_dim, hidden_dim)
        self.dna_encoder = SymmetricDNAEncoder(hidden_dim=hidden_dim)
        self.comparative = _mlp(hidden_dim * 4, hidden_dim, hidden_dim)
        self.film = nn.Linear(hidden_dim, hidden_dim * 2)
        self.head = nn.Linear(hidden_dim, 1)

    def forward(
        self,
        protein_embedding: torch.Tensor,
        target_window_features_batch: torch.Tensor,
        candidate_features: torch.Tensor,
    ) -> torch.Tensor:
        protein = self.protein_projection(protein_embedding)
        target = self.dna_encoder(target_window_features_batch).mean(dim=1)
        candidate = self.dna_encoder(candidate_features)
        comparative = torch.cat(
            [candidate, target, candidate - target, candidate * target],
            dim=-1,
        )
        hidden = self.comparative(comparative)
        gamma, beta = self.film(protein).chunk(2, dim=-1)
        hidden = (1.0 + torch.tanh(gamma)) * hidden + beta
        return self.head(hidden).squeeze(-1)


def model_parameter_counts(model: nn.Module) -> tuple[int, int]:
    trainable = sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)
    frozen = sum(parameter.numel() for parameter in model.parameters() if not parameter.requires_grad)
    return trainable, frozen


def score_rc(
    model: nn.Module,
    model_name: str,
    candidate_sequences: Sequence[str],
    *,
    protein_embedding: torch.Tensor | None = None,
    target: str | None = None,
    device: str = "cpu",
) -> np.ndarray:
    """Score sequences using the same RC-symmetric representation for every model."""
    candidates = batch_rc_symmetric_one_hot(candidate_sequences).to(device)
    target_features = None
    if target is not None:
        target_features = target_window_features(target).to(device)
    model = model.to(device)
    model.eval()
    with torch.no_grad():
        if model_name == CandidateDNAOnly.model_name:
            output = model(candidates)
        elif model_name == ProteinCandidate.model_name:
            if protein_embedding is None:
                raise ValueError("M1 requires a frozen protein embedding")
            protein = protein_embedding.to(device).reshape(1, -1).expand(len(candidate_sequences), -1)
            output = model(protein, candidates)
        elif model_name == TargetCandidateOnly.model_name:
            if target_features is None:
                raise ValueError("M2 requires a target")
            target_batch = target_features.unsqueeze(0).expand(len(candidate_sequences), -1, -1)
            output = model(target_batch, candidates)
        elif model_name == ProteinTargetCandidate.model_name:
            if protein_embedding is None or target_features is None:
                raise ValueError("M3 requires a frozen protein embedding and target")
            protein = protein_embedding.to(device).reshape(1, -1).expand(len(candidate_sequences), -1)
            target_batch = target_features.unsqueeze(0).expand(len(candidate_sequences), -1, -1)
            output = model(protein, target_batch, candidates)
        else:
            raise ValueError(f"Unknown model name: {model_name}")
    return output.detach().cpu().numpy()


def _best_window_score(target: str, candidate: str, scorer) -> float:
    candidate = validate_dna_sequence(candidate, min_length=KMER_LENGTH)
    orientations = (candidate, reverse_complement(candidate))
    return max(float(scorer(window, oriented)) for window in target_windows(target) for oriented in orientations)


def target_hamming_control(target: str, candidate: str) -> float:
    return _best_window_score(
        target,
        candidate,
        lambda window, oriented: 1.0
        - sum(a != b for a, b in zip(window, oriented)) / len(window),
    )


def target_edit_control(target: str, candidate: str) -> float:
    return _best_window_score(
        target,
        candidate,
        lambda window, oriented: 1.0
        - edit_distance(window, oriented) / max(len(window), len(oriented)),
    )


def target_kmer_overlap_control(target: str, candidate: str, k: int = 3) -> float:
    return _best_window_score(
        target,
        candidate,
        lambda window, oriented: kmer_jaccard(window, oriented, k=k, rc_aware=True),
    )


def toy_nonseparable_score(protein_scalar: float, target_scalar: float, candidate_scalar: float) -> float:
    """Toy counterexample used only to test expressible target-dependent ranking reversal."""
    return protein_scalar * target_scalar * candidate_scalar
