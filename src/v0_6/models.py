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
DNA_FEATURE_DIM = KMER_LENGTH * 4
TARGET_FEATURE_DIM = DNA_FEATURE_DIM + 1


def validate_dna_sequence(sequence: str, *, min_length: int = 1) -> str:
    normalized = normalize_sequence(sequence)
    if len(normalized) < min_length or any(base not in DNA_INDEX for base in normalized):
        raise ValueError(f"Expected ACGT DNA of length >= {min_length}, got {sequence!r}")
    return normalized


def canonical_sequence(sequence: str) -> str:
    normalized = validate_dna_sequence(sequence)
    return canonical_rc(normalized)


def _one_hot(sequence: str) -> np.ndarray:
    sequence = validate_dna_sequence(sequence)
    result = np.zeros((len(sequence), len(DNA_ALPHABET)), dtype=np.float32)
    for index, base in enumerate(sequence):
        result[index, DNA_INDEX[base]] = 1.0
    return result.reshape(-1)


def canonical_one_hot(sequence: str) -> np.ndarray:
    """Injective representation of an RC class using its canonical strand."""
    return _one_hot(canonical_sequence(sequence)).astype(np.float32)


def batch_canonical_one_hot(sequences: Sequence[str]) -> torch.Tensor:
    if not sequences:
        raise ValueError("At least one DNA sequence is required")
    lengths = {len(validate_dna_sequence(sequence)) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("A batch must contain equal-length DNA sequences")
    return torch.from_numpy(np.vstack([canonical_one_hot(sequence) for sequence in sequences]))


def ordered_target_windows(
    target: str,
    k: int = KMER_LENGTH,
) -> tuple[tuple[str, int, float], ...]:
    """Return every target window in canonical-strand order with its register.

    Canonicalizing the complete target makes the representation invariant to a
    supplied target orientation. Windows are deliberately not canonicalized,
    deduplicated, or sorted independently, so order and repeated motifs stay
    observable by the window encoder.
    """
    target = canonical_sequence(validate_dna_sequence(target, min_length=k))
    n_windows = len(target) - k + 1
    if n_windows < 1:
        raise ValueError(f"Target must contain at least one {k}-mer window")
    denominator = max(1, n_windows - 1)
    return tuple(
        (
            target[index : index + k],
            index,
            (2.0 * index / denominator) - 1.0 if n_windows > 1 else 0.0,
        )
        for index in range(n_windows)
    )


def target_window_features_v06(target: str, k: int = KMER_LENGTH) -> torch.Tensor:
    """Encode ordered windows as 28 DNA features plus one register feature."""
    rows = []
    for window, _index, position in ordered_target_windows(target, k=k):
        # The full target is canonicalized once; each ordered window retains
        # its strand/orientation relative to that canonical target.
        rows.append(np.concatenate([_one_hot(window), np.asarray([position], dtype=np.float32)]))
    return torch.from_numpy(np.stack(rows).astype(np.float32))


def _pool_target_windows(values: torch.Tensor, positions: torch.Tensor, mode: str) -> torch.Tensor:
    """Small deterministic pooling operators used by the audit."""
    if values.ndim != 2:
        raise ValueError(f"Expected [windows,hidden] values, got {tuple(values.shape)}")
    if mode == "mean":
        return values.mean(dim=0)
    if mode == "max":
        return values.max(dim=0).values
    if mode in {"logsumexp", "lse"}:
        return torch.logsumexp(values, dim=0) - np.log(values.shape[0])
    if mode == "softmax":
        weights = torch.softmax(values.mean(dim=1), dim=0).unsqueeze(-1)
        return (values * weights).sum(dim=0)
    if mode in {"position_aware", "position_weighted_mean"}:
        if positions.ndim != 1 or positions.shape[0] != values.shape[0]:
            raise ValueError("Positions must align with target windows")
        weights = (1.0 + positions.abs()).unsqueeze(-1)
        return (values * weights).sum(dim=0) / weights.sum()
    raise ValueError(f"Unknown target pooling mode: {mode}")


class CandidateDNAEncoderV06(nn.Module):
    def __init__(self, hidden_dim: int = 32) -> None:
        super().__init__()
        self.projection = nn.Sequential(nn.Linear(DNA_FEATURE_DIM, hidden_dim), nn.Tanh())
        self.hidden_dim = hidden_dim

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        if features.ndim != 2 or features.shape[-1] != DNA_FEATURE_DIM:
            raise ValueError(f"Expected [batch,{DNA_FEATURE_DIM}], got {tuple(features.shape)}")
        return self.projection(features)


class OrderedTargetEncoderV06(nn.Module):
    def __init__(self, hidden_dim: int = 32) -> None:
        super().__init__()
        self.projection = nn.Sequential(nn.Linear(TARGET_FEATURE_DIM, hidden_dim), nn.Tanh())
        self.hidden_dim = hidden_dim

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        if features.ndim != 3 or features.shape[-1] != TARGET_FEATURE_DIM:
            raise ValueError(
                f"Expected [batch,windows,{TARGET_FEATURE_DIM}], got {tuple(features.shape)}"
            )
        return self.projection(features)


class ProteinProjectionV06(nn.Module):
    def __init__(self, protein_dim: int, hidden_dim: int) -> None:
        super().__init__()
        self.projection = nn.Sequential(nn.Linear(protein_dim, hidden_dim), nn.Tanh())

    def forward(self, protein_embedding: torch.Tensor) -> torch.Tensor:
        return self.projection(protein_embedding)


def _mlp(input_dim: int, hidden_dim: int, output_dim: int = 1) -> nn.Sequential:
    return nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, output_dim))


class CandidateDNAOnlyV06(nn.Module):
    model_name = "M0_CandidateDNAOnly_v0_6"
    inputs = "D"

    def __init__(self, hidden_dim: int = 32) -> None:
        super().__init__()
        self.dna_encoder = CandidateDNAEncoderV06(hidden_dim)
        self.head = _mlp(hidden_dim, hidden_dim)

    def forward(self, candidate_features: torch.Tensor) -> torch.Tensor:
        return self.head(self.dna_encoder(candidate_features)).squeeze(-1)


class ProteinCandidateV06(nn.Module):
    model_name = "M1_ProteinCandidate_v0_6"
    inputs = "P,D"

    def __init__(self, protein_dim: int, hidden_dim: int = 32) -> None:
        super().__init__()
        self.protein_projection = ProteinProjectionV06(protein_dim, hidden_dim)
        self.dna_encoder = CandidateDNAEncoderV06(hidden_dim)
        self.head = _mlp(hidden_dim * 3, hidden_dim)

    def forward(self, protein_embedding: torch.Tensor, candidate_features: torch.Tensor) -> torch.Tensor:
        protein = self.protein_projection(protein_embedding)
        candidate = self.dna_encoder(candidate_features)
        return self.head(torch.cat([protein, candidate, protein * candidate], dim=-1)).squeeze(-1)


class _TargetComparativeBaseV06(nn.Module):
    def __init__(self, hidden_dim: int = 32, target_pooling: str = "mean") -> None:
        super().__init__()
        self.dna_encoder = CandidateDNAEncoderV06(hidden_dim)
        self.target_encoder = OrderedTargetEncoderV06(hidden_dim)
        self.comparative = _mlp(hidden_dim * 4, hidden_dim, hidden_dim)
        self.target_pooling = target_pooling

    def _comparative_hidden(
        self,
        target_window_features_batch: torch.Tensor,
        candidate_features: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        candidate = self.dna_encoder(candidate_features)
        target_windows = self.target_encoder(target_window_features_batch)
        if target_windows.shape[0] != candidate.shape[0]:
            raise ValueError("Target and candidate batches must have matching batch dimensions")
        comparisons = torch.cat(
            [
                candidate.unsqueeze(1).expand_as(target_windows),
                target_windows,
                candidate.unsqueeze(1) - target_windows,
                candidate.unsqueeze(1) * target_windows,
            ],
            dim=-1,
        )
        hidden = self.comparative(comparisons)
        positions = target_window_features_batch[..., -1]
        return hidden, positions

    def pool(self, values: torch.Tensor, positions: torch.Tensor) -> torch.Tensor:
        return torch.stack(
            [_pool_target_windows(value, position, self.target_pooling) for value, position in zip(values, positions)]
        )


class TargetCandidateOnlyV06(_TargetComparativeBaseV06):
    model_name = "M2_TargetCandidateOnly_v0_6"
    inputs = "T,D"

    def __init__(self, hidden_dim: int = 32, target_pooling: str = "mean") -> None:
        super().__init__(hidden_dim, target_pooling)
        self.head = nn.Linear(hidden_dim, 1)

    def forward(self, target_window_features_batch: torch.Tensor, candidate_features: torch.Tensor) -> torch.Tensor:
        hidden, positions = self._comparative_hidden(target_window_features_batch, candidate_features)
        return self.head(self.pool(hidden, positions)).squeeze(-1)


class ProteinTargetCandidateV06(_TargetComparativeBaseV06):
    model_name = "M3_ProteinTargetCandidate_v0_6"
    inputs = "P,T,D"

    def __init__(self, protein_dim: int, hidden_dim: int = 32, target_pooling: str = "mean") -> None:
        super().__init__(hidden_dim, target_pooling)
        self.protein_projection = ProteinProjectionV06(protein_dim, hidden_dim)
        self.film = nn.Linear(hidden_dim, hidden_dim * 2)
        self.head = nn.Linear(hidden_dim, 1)

    def forward(
        self,
        protein_embedding: torch.Tensor,
        target_window_features_batch: torch.Tensor,
        candidate_features: torch.Tensor,
    ) -> torch.Tensor:
        hidden, positions = self._comparative_hidden(target_window_features_batch, candidate_features)
        protein = self.protein_projection(protein_embedding).unsqueeze(1)
        gamma, beta = self.film(protein).chunk(2, dim=-1)
        hidden = (1.0 + torch.tanh(gamma)) * hidden + beta
        return self.head(self.pool(hidden, positions)).squeeze(-1)


MODEL_ORDER_V06 = ("M0", "M1", "M1c", "M2", "M3")
MODEL_FULL_NAMES_V06 = {
    "M0": CandidateDNAOnlyV06.model_name,
    "M1": ProteinCandidateV06.model_name,
    "M1c": ProteinCandidateV06.model_name,
    "M2": TargetCandidateOnlyV06.model_name,
    "M3": ProteinTargetCandidateV06.model_name,
}


def model_parameter_counts_v06(model: nn.Module) -> tuple[int, int]:
    trainable = sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)
    frozen = sum(parameter.numel() for parameter in model.parameters() if not parameter.requires_grad)
    return trainable, frozen


def score_rc_v06(
    model: nn.Module,
    model_name: str,
    candidate_sequences: Sequence[str],
    *,
    protein_embedding: torch.Tensor | None = None,
    target: str | None = None,
    target_pooling: str = "mean",
    device: str = "cpu",
) -> np.ndarray:
    candidates = batch_canonical_one_hot(candidate_sequences).to(device)
    model = model.to(device)
    model.eval()
    target_features = None
    if target is not None:
        target_features = target_window_features_v06(target).to(device)
    with torch.no_grad():
        if model_name == CandidateDNAOnlyV06.model_name:
            output = model(candidates)
        elif model_name == ProteinCandidateV06.model_name:
            if protein_embedding is None:
                raise ValueError("M1 requires a frozen protein embedding")
            protein = protein_embedding.to(device).reshape(1, -1).expand(len(candidate_sequences), -1)
            output = model(protein, candidates)
        elif model_name == TargetCandidateOnlyV06.model_name:
            if target_features is None:
                raise ValueError("M2 requires a target")
            target_batch = target_features.unsqueeze(0).expand(len(candidate_sequences), -1, -1)
            output = model(target_batch, candidates)
        elif model_name == ProteinTargetCandidateV06.model_name:
            if protein_embedding is None or target_features is None:
                raise ValueError("M3 requires a frozen protein embedding and target")
            protein = protein_embedding.to(device).reshape(1, -1).expand(len(candidate_sequences), -1)
            target_batch = target_features.unsqueeze(0).expand(len(candidate_sequences), -1, -1)
            output = model(protein, target_batch, candidates)
        else:
            raise ValueError(f"Unknown v0.6 model name: {model_name}")
    return output.detach().cpu().numpy()


def _best_window_score(target: str, candidate: str, scorer) -> float:
    candidate = validate_dna_sequence(candidate, min_length=KMER_LENGTH)
    orientations = (candidate, reverse_complement(candidate))
    return max(
        float(scorer(window, oriented))
        for window, _index, _position in ordered_target_windows(target)
        for oriented in orientations
    )


def target_hamming_control_v06(target: str, candidate: str) -> float:
    return _best_window_score(
        target,
        candidate,
        lambda window, oriented: 1.0 - sum(a != b for a, b in zip(window, oriented)) / len(window),
    )


def target_edit_control_v06(target: str, candidate: str) -> float:
    return _best_window_score(
        target,
        candidate,
        lambda window, oriented: 1.0 - edit_distance(window, oriented) / max(len(window), len(oriented)),
    )


def target_kmer_overlap_control_v06(target: str, candidate: str, k: int = 3) -> float:
    return _best_window_score(
        target,
        candidate,
        lambda window, oriented: kmer_jaccard(window, oriented, k=k, rc_aware=True),
    )
