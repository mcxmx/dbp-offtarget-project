from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import torch

from src.sequence_equivalence import reverse_complement
from src.v0_5_local_models import (
    LocalProteinCandidate,
    LocalProteinCandidateCapacityMatched,
    LocalProteinTargetCandidate,
    batch_rc_oriented_dna,
    local_model_parameter_counts,
    score_local_rc,
    target_window_oriented_features,
)
from src.v0_5_models import toy_nonseparable_score


ROOT = Path(__file__).resolve().parents[1]
PROTEIN_DIM = 480


def _inputs():
    torch.manual_seed(7)
    protein = torch.randn(9, PROTEIN_DIM)
    candidates = batch_rc_oriented_dna(["ACGTCAG", "TTTTTTT"])
    target_a = target_window_oriented_features("TAGCAGGATGTGT")
    target_b = target_window_oriented_features("GCAGATCTGCACATC")
    return protein, candidates, target_a, target_b


def test_residue_attention_cache_shape_contract():
    protein, candidates, _, _ = _inputs()
    model = LocalProteinCandidate(PROTEIN_DIM, hidden_dim=8)
    score, attention = model.forward_with_attention(protein, candidates)
    assert score.shape == (2,)
    assert attention.shape == (2, 7, 9)


def test_candidate_dependent_residue_attention():
    protein, candidates, _, _ = _inputs()
    model = LocalProteinCandidate(PROTEIN_DIM, hidden_dim=8)
    _, first_attention = model.forward_with_attention(protein, candidates[:1])
    _, second_attention = model.forward_with_attention(protein, candidates[1:2])
    assert not torch.allclose(first_attention, second_attention)


def test_local_model_forward_contracts():
    protein, candidates, target_a, _ = _inputs()
    l1 = LocalProteinCandidate(PROTEIN_DIM, hidden_dim=8)
    l1c = LocalProteinCandidateCapacityMatched(PROTEIN_DIM, hidden_dim=8, head_hidden_dim=16)
    l2 = LocalProteinTargetCandidate(PROTEIN_DIM, hidden_dim=8)
    assert l1(protein, candidates).shape == (2,)
    assert l1c(protein, candidates).shape == (2,)
    assert l2(protein, target_a, candidates).shape == (2,)


def test_l1_prediction_changes_with_protein():
    _, candidates, _, _ = _inputs()
    model = LocalProteinCandidate(PROTEIN_DIM, hidden_dim=8)
    first = model(torch.zeros(9, PROTEIN_DIM), candidates)
    second = model(torch.ones(9, PROTEIN_DIM), candidates)
    assert not torch.allclose(first, second)


def test_l2_prediction_changes_with_protein_and_target():
    protein, candidates, target_a, target_b = _inputs()
    model = LocalProteinTargetCandidate(PROTEIN_DIM, hidden_dim=8)
    first = model(protein, target_a, candidates)
    target_changed = model(protein, target_b, candidates)
    protein_changed = model(torch.zeros_like(protein), target_a, candidates)
    assert not torch.allclose(first, target_changed)
    assert not torch.allclose(first, protein_changed)


def test_l2_has_target_dependent_nonseparable_expression():
    # The same toy counterexample used by the frozen P,T,D contract demonstrates
    # that target-dependent ranking reversal is expressible in the hypothesis.
    assert toy_nonseparable_score(1.0, 1.0, 1.0) > toy_nonseparable_score(1.0, 1.0, -1.0)
    assert toy_nonseparable_score(1.0, -1.0, 1.0) < toy_nonseparable_score(1.0, -1.0, -1.0)


def test_local_candidate_rc_invariance():
    protein, _, _, _ = _inputs()
    model = LocalProteinCandidate(PROTEIN_DIM, hidden_dim=8)
    sequence = "ACGTCAG"
    first = score_local_rc(model, model.model_name, protein, [sequence])
    second = score_local_rc(model, model.model_name, protein, [reverse_complement(sequence)])
    assert np.allclose(first, second, atol=1e-6)


def test_local_target_orientation_invariance():
    protein, _, _, _ = _inputs()
    model = LocalProteinTargetCandidate(PROTEIN_DIM, hidden_dim=8)
    sequence = "ACGTCAG"
    target = "TAGCAGGATGTGT"
    first = score_local_rc(model, model.model_name, protein, [sequence], target=target)
    second = score_local_rc(
        model,
        model.model_name,
        protein,
        [sequence],
        target=reverse_complement(target),
    )
    assert np.allclose(first, second, atol=1e-6)


def test_local_parameter_counts_are_small_and_capacity_matched():
    l1 = local_model_parameter_counts(LocalProteinCandidate(PROTEIN_DIM, hidden_dim=24))[0]
    l1c = local_model_parameter_counts(
        LocalProteinCandidateCapacityMatched(PROTEIN_DIM, hidden_dim=24, head_hidden_dim=64)
    )[0]
    l2 = local_model_parameter_counts(LocalProteinTargetCandidate(PROTEIN_DIM, hidden_dim=24))[0]
    assert max(l1, l1c, l2) < 100_000
    assert l2 / l1c < 2.0
    assert l2 / l1c > 0.5


def test_input_signatures_keep_model_blindness_contracts():
    assert list(inspect.signature(LocalProteinCandidate.forward).parameters) == [
        "self",
        "protein_residues",
        "candidate_features",
    ]
    assert list(inspect.signature(LocalProteinTargetCandidate.forward).parameters) == [
        "self",
        "protein_residues",
        "target_features",
        "candidate_features",
    ]
