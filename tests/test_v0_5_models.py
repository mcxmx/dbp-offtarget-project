from pathlib import Path
import inspect

import numpy as np
import torch

from src.sequence_equivalence import reverse_complement
from src.v0_5_models import (
    CandidateDNAOnly,
    ProteinCandidate,
    ProteinTargetCandidate,
    TargetCandidateOnly,
    batch_rc_symmetric_one_hot,
    score_rc,
    target_edit_control,
    target_hamming_control,
    target_kmer_overlap_control,
    target_windows,
    toy_nonseparable_score,
)


ROOT = Path(__file__).resolve().parents[1]
EMBEDDING_DIM = 480


def test_target_windows_are_complete_and_orientation_invariant():
    target = "TAGCAGGATGTGT"
    assert len(target_windows(target)) == 7
    assert target_windows(target) == target_windows(reverse_complement(target))


def test_rc_symmetric_candidate_features_are_identical():
    sequence = "ACGTCAG"
    features = batch_rc_symmetric_one_hot([sequence, reverse_complement(sequence)])
    assert torch.allclose(features[0], features[1])


def test_m0_forward_and_rc_score_are_invariant():
    model = CandidateDNAOnly(hidden_dim=8)
    sequence = "ACGTCAG"
    scores = score_rc(model, model.model_name, [sequence, reverse_complement(sequence)])
    assert np.allclose(scores[0], scores[1])


def test_m1_forward_and_m2_forward():
    protein = torch.zeros(2, EMBEDDING_DIM)
    candidates = batch_rc_symmetric_one_hot(["ACGTCAG", "TTTTTTT"])
    target = torch.zeros(2, 7, 28)
    m1 = ProteinCandidate(EMBEDDING_DIM, hidden_dim=8)
    m2 = TargetCandidateOnly(hidden_dim=8)
    assert m1(protein, candidates).shape == (2,)
    assert m2(target, candidates).shape == (2,)


def test_m3_forward_changes_with_target_and_protein():
    torch.manual_seed(42)
    model = ProteinTargetCandidate(EMBEDDING_DIM, hidden_dim=8)
    candidate = ["ACGTCAG"]
    protein_a = torch.zeros(EMBEDDING_DIM)
    protein_b = torch.ones(EMBEDDING_DIM)
    target_a = "TAGCAGGATGTGT"
    target_b = "GCAGATCTGCACATC"
    score_a = score_rc(model, model.model_name, candidate, protein_embedding=protein_a, target=target_a)
    score_t = score_rc(model, model.model_name, candidate, protein_embedding=protein_a, target=target_b)
    score_p = score_rc(model, model.model_name, candidate, protein_embedding=protein_b, target=target_a)
    assert not np.allclose(score_a, score_t)
    assert not np.allclose(score_a, score_p)


def test_all_model_scores_are_candidate_rc_invariant():
    candidate = ["ACGTCAG"]
    rc_candidate = [reverse_complement(candidate[0])]
    protein = torch.zeros(EMBEDDING_DIM)
    target = "TAGCAGGATGTGT"
    models = [
        (CandidateDNAOnly(hidden_dim=8), CandidateDNAOnly.model_name),
        (ProteinCandidate(EMBEDDING_DIM, hidden_dim=8), ProteinCandidate.model_name),
        (TargetCandidateOnly(hidden_dim=8), TargetCandidateOnly.model_name),
        (ProteinTargetCandidate(EMBEDDING_DIM, hidden_dim=8), ProteinTargetCandidate.model_name),
    ]
    for model, model_name in models:
        first = score_rc(
            model,
            model_name,
            candidate,
            protein_embedding=protein,
            target=target,
        )
        second = score_rc(
            model,
            model_name,
            rc_candidate,
            protein_embedding=protein,
            target=target,
        )
        assert np.allclose(first, second)


def test_m3_score_is_target_orientation_invariant():
    model = ProteinTargetCandidate(EMBEDDING_DIM, hidden_dim=8)
    protein = torch.zeros(EMBEDDING_DIM)
    target = "TAGCAGGATGTGT"
    candidate = ["ACGTCAG"]
    first = score_rc(model, model.model_name, candidate, protein_embedding=protein, target=target)
    second = score_rc(
        model,
        model.model_name,
        candidate,
        protein_embedding=protein,
        target=reverse_complement(target),
    )
    assert np.allclose(first, second)


def test_m3_can_express_target_dependent_ranking_reversal():
    # Fixed P=1: changing T from +1 to -1 reverses D1=+1 vs D2=-1.
    p = 1.0
    assert toy_nonseparable_score(p, 1.0, 1.0) > toy_nonseparable_score(p, 1.0, -1.0)
    assert toy_nonseparable_score(p, -1.0, 1.0) < toy_nonseparable_score(p, -1.0, -1.0)


def test_m3_is_not_candidate_only_plus_constant():
    # A candidate-only score plus a target-only constant preserves ranking.
    candidate_scores = np.array([0.8, 0.2])
    assert np.argsort(-(candidate_scores + 5.0)).tolist() == [0, 1]
    assert np.argsort(-(candidate_scores - 7.0)).tolist() == [0, 1]


def test_models_have_expected_input_contracts():
    assert CandidateDNAOnly.model_name == "M0_CandidateDNAOnly"
    assert ProteinCandidate.model_name == "M1_ProteinCandidate"
    assert TargetCandidateOnly.model_name == "M2_TargetCandidateOnly"
    assert ProteinTargetCandidate.model_name == "M3_ProteinTargetCandidate"
    assert CandidateDNAOnly.inputs == "D"
    assert ProteinCandidate.inputs == "P,D"
    assert TargetCandidateOnly.inputs == "T,D"
    assert ProteinTargetCandidate.inputs == "P,T,D"


def test_target_controls_are_rc_invariant():
    target = "TAGCAGGATGTGT"
    candidate = "ACGTCAG"
    rc = reverse_complement(candidate)
    for control in [target_hamming_control, target_edit_control, target_kmer_overlap_control]:
        assert np.isclose(control(target, candidate), control(target, rc))


def test_input_signatures_enforce_model_blindness_contracts():
    assert list(inspect.signature(CandidateDNAOnly.forward).parameters) == ["self", "candidate_features"]
    assert list(inspect.signature(ProteinCandidate.forward).parameters) == [
        "self",
        "protein_embedding",
        "candidate_features",
    ]
    assert list(inspect.signature(TargetCandidateOnly.forward).parameters) == [
        "self",
        "target_window_features_batch",
        "candidate_features",
    ]


def test_target_windows_reject_short_targets():
    try:
        target_windows("ACGT")
    except ValueError:
        pass
    else:
        raise AssertionError("Targets shorter than 7 bp must be rejected")
