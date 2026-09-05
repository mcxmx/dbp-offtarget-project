from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import torch
from torch import nn

from src.sequence_equivalence import reverse_complement
from src.v0_5_models import DNA_ALPHABET, DNA_INDEX, KMER_LENGTH, target_windows, validate_dna_sequence


def oriented_one_hot(sequence: str) -> np.ndarray:
    sequence = validate_dna_sequence(sequence, min_length=KMER_LENGTH)
    if len(sequence) != KMER_LENGTH:
        raise ValueError(f"Expected a {KMER_LENGTH}-mer, got {sequence!r}")
    features = np.zeros((KMER_LENGTH, len(DNA_ALPHABET)), dtype=np.float32)
    for index, base in enumerate(sequence):
        features[index, DNA_INDEX[base]] = 1.0
    return features


def batch_rc_oriented_dna(sequences: Sequence[str]) -> torch.Tensor:
    if not sequences:
        raise ValueError("At least one DNA sequence is required")
    return torch.from_numpy(
        np.stack(
            [
                np.stack([oriented_one_hot(sequence), oriented_one_hot(reverse_complement(sequence))])
                for sequence in sequences
            ]
        )
    )


def target_window_oriented_features(target: str, k: int = KMER_LENGTH) -> torch.Tensor:
    windows = target_windows(target, k=k)
    return batch_rc_oriented_dna(windows)


def _as_oriented_features(features: torch.Tensor) -> torch.Tensor:
    if features.ndim == 3:
        features = features.unsqueeze(1)
    if features.ndim != 4 or tuple(features.shape[-2:]) != (KMER_LENGTH, len(DNA_ALPHABET)):
        raise ValueError(
            "Expected DNA features [batch,orientations,7,4] or [batch,7,4], "
            f"got {tuple(features.shape)}"
        )
    if features.shape[1] != 2:
        raise ValueError(f"Expected two orientations, got {features.shape[1]}")
    return features.float()


class ResidueDNAAttention(nn.Module):
    """One small candidate-dependent attention block over frozen protein residues."""

    def __init__(self, protein_dim: int = 480, hidden_dim: int = 24) -> None:
        super().__init__()
        self.protein_dim = protein_dim
        self.hidden_dim = hidden_dim
        self.protein_projection = nn.Linear(protein_dim, hidden_dim)
        self.dna_projection = nn.Linear(len(DNA_ALPHABET), hidden_dim)
        self.position_embedding = nn.Parameter(torch.zeros(KMER_LENGTH, hidden_dim))
        nn.init.normal_(self.position_embedding, mean=0.0, std=0.02)
        self.interaction_mlp = nn.Sequential(
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.Tanh(),
        )

    def _expand_protein(self, protein_residues: torch.Tensor, batch_size: int) -> torch.Tensor:
        if protein_residues.ndim == 2:
            return protein_residues.unsqueeze(0).expand(batch_size, -1, -1).float()
        if protein_residues.ndim == 3:
            if protein_residues.shape[0] == 1:
                return protein_residues.expand(batch_size, -1, -1).float()
            if protein_residues.shape[0] == batch_size:
                return protein_residues.float()
        raise ValueError(
            "Expected protein residue embeddings [length,dim] or [batch,length,dim], "
            f"got {tuple(protein_residues.shape)}"
        )

    def _forward_oriented(
        self,
        protein_residues: torch.Tensor,
        dna_features: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        protein = self.protein_projection(protein_residues)
        dna = self.dna_projection(dna_features) + self.position_embedding.unsqueeze(0)
        logits = torch.einsum("bih,brh->bir", dna, protein) / np.sqrt(self.hidden_dim)
        attention = torch.softmax(logits, dim=-1)
        context = torch.einsum("bir,brh->bih", attention, protein)
        comparison = torch.cat([dna, context, dna - context, dna * context], dim=-1)
        position_hidden = self.interaction_mlp(comparison)
        return position_hidden.mean(dim=1), attention

    def forward_rc(
        self,
        protein_residues: torch.Tensor,
        dna_features: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        oriented = _as_oriented_features(dna_features)
        batch_size, orientations = oriented.shape[:2]
        protein = self._expand_protein(protein_residues, batch_size)
        flat_dna = oriented.reshape(batch_size * orientations, KMER_LENGTH, len(DNA_ALPHABET))
        flat_protein = (
            protein.unsqueeze(1)
            .expand(batch_size, orientations, protein.shape[1], protein.shape[2])
            .reshape(batch_size * orientations, protein.shape[1], protein.shape[2])
        )
        vectors, attention = self._forward_oriented(flat_protein, flat_dna)
        vectors = vectors.reshape(batch_size, orientations, self.hidden_dim).mean(dim=1)
        attention = attention.reshape(batch_size, orientations, KMER_LENGTH, protein.shape[1]).mean(dim=1)
        return vectors, attention

    def forward(
        self,
        protein_residues: torch.Tensor,
        dna_features: torch.Tensor,
    ) -> torch.Tensor:
        return self.forward_rc(protein_residues, dna_features)[0]


def _small_mlp(input_dim: int, hidden_dim: int, output_dim: int = 1) -> nn.Sequential:
    return nn.Sequential(
        nn.Linear(input_dim, hidden_dim),
        nn.Tanh(),
        nn.Linear(hidden_dim, output_dim),
    )


class LocalProteinCandidate(nn.Module):
    model_name = "L1_LocalProteinCandidate"
    inputs = "P,D"

    def __init__(self, protein_dim: int = 480, hidden_dim: int = 24, head_hidden_dim: int | None = None) -> None:
        super().__init__()
        head_hidden_dim = hidden_dim if head_hidden_dim is None else head_hidden_dim
        self.interaction = ResidueDNAAttention(protein_dim=protein_dim, hidden_dim=hidden_dim)
        self.head = _small_mlp(hidden_dim, head_hidden_dim)

    def forward_with_attention(
        self,
        protein_residues: torch.Tensor,
        candidate_features: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        interaction, attention = self.interaction.forward_rc(protein_residues, candidate_features)
        return self.head(interaction).squeeze(-1), attention

    def forward(self, protein_residues: torch.Tensor, candidate_features: torch.Tensor) -> torch.Tensor:
        return self.forward_with_attention(protein_residues, candidate_features)[0]


class LocalProteinCandidateCapacityMatched(LocalProteinCandidate):
    model_name = "L1c_LocalProteinCandidateCapacityMatched"

    def __init__(self, protein_dim: int = 480, hidden_dim: int = 24, head_hidden_dim: int = 64) -> None:
        super().__init__(protein_dim=protein_dim, hidden_dim=hidden_dim, head_hidden_dim=head_hidden_dim)


class LocalProteinTargetCandidate(nn.Module):
    model_name = "L2_LocalProteinTargetCandidate"
    inputs = "P,T,D"

    def __init__(self, protein_dim: int = 480, hidden_dim: int = 24) -> None:
        super().__init__()
        self.hidden_dim = hidden_dim
        self.interaction = ResidueDNAAttention(protein_dim=protein_dim, hidden_dim=hidden_dim)
        self.comparative = _small_mlp(hidden_dim * 4, hidden_dim, hidden_dim)
        self.head = nn.Linear(hidden_dim, 1)

    def forward_with_attention(
        self,
        protein_residues: torch.Tensor,
        target_features: torch.Tensor,
        candidate_features: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        candidate, candidate_attention = self.interaction.forward_rc(protein_residues, candidate_features)
        target, target_attention = self.interaction.forward_rc(protein_residues, target_features)
        if target.ndim != 2:
            raise ValueError("Target features must represent a batch of target windows")
        candidate_expanded = candidate.unsqueeze(1).expand(-1, target.shape[0], -1)
        target_expanded = target.unsqueeze(0).expand(candidate.shape[0], -1, -1)
        comparative = torch.cat(
            [
                candidate_expanded,
                target_expanded,
                candidate_expanded - target_expanded,
                candidate_expanded * target_expanded,
            ],
            dim=-1,
        )
        window_hidden = self.comparative(comparative)
        window_scores = self.head(window_hidden).squeeze(-1)
        return window_scores.mean(dim=1), candidate_attention, target_attention

    def forward(
        self,
        protein_residues: torch.Tensor,
        target_features: torch.Tensor,
        candidate_features: torch.Tensor,
    ) -> torch.Tensor:
        return self.forward_with_attention(protein_residues, target_features, candidate_features)[0]


def score_local_rc(
    model: nn.Module,
    model_name: str,
    protein_residues: torch.Tensor,
    candidate_sequences: Sequence[str],
    *,
    target: str | None = None,
    device: str = "cpu",
) -> np.ndarray:
    candidates = batch_rc_oriented_dna(candidate_sequences).to(device)
    protein_residues = protein_residues.to(device)
    model = model.to(device)
    model.eval()
    with torch.no_grad():
        if model_name in {LocalProteinCandidate.model_name, LocalProteinCandidateCapacityMatched.model_name}:
            scores = model(protein_residues, candidates)
        elif model_name == LocalProteinTargetCandidate.model_name:
            if target is None:
                raise ValueError("L2 requires a primary target")
            target_features = target_window_oriented_features(target).to(device)
            scores = model(protein_residues, target_features, candidates)
        else:
            raise ValueError(f"Unknown local model name: {model_name}")
    return scores.detach().cpu().numpy()


def local_model_parameter_counts(model: nn.Module) -> tuple[int, int]:
    trainable = sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)
    frozen = sum(parameter.numel() for parameter in model.parameters() if not parameter.requires_grad)
    return trainable, frozen
