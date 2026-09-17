from itertools import product
from pathlib import Path

import numpy as np
import torch

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.v0_5_models import rc_symmetric_one_hot
from src.v0_6.audit import all_kmers, check_score_rc_invariance, rc_representation_audit, synthetic_conditionality_benchmark
from src.v0_6.gates import assert_pre_replay_gates
from src.v0_6.models import (
    CandidateDNAOnlyV06,
    ProteinCandidateV06,
    ProteinTargetCandidateV06,
    TargetCandidateOnlyV06,
    batch_canonical_one_hot,
    canonical_one_hot,
    ordered_target_windows,
    score_rc_v06,
    target_window_features_v06,
)


def test_exhaustive_rc_class_count_and_v06_injectivity():
    audit = rc_representation_audit().iloc[0]
    assert audit["n_oriented_7mers"] == 16384
    assert audit["n_canonical_rc_classes"] == 8192
    assert audit["legacy_representation_count"] == 2000
    assert audit["legacy_cross_rc_class_representation_collisions"] == 1872
    assert audit["v06_representation_count"] == 8192
    assert audit["v06_cross_rc_class_representation_collisions"] == 0
    assert bool(audit["v06_rc_class_injective"])


def test_v06_representations_equal_iff_rc_equivalent():
    sequences = all_kmers()
    representations = {}
    for sequence in sequences:
        representations.setdefault(tuple(canonical_one_hot(sequence)), set()).add(canonical_rc(sequence))
    assert len(representations) == 8192
    assert all(len(classes) == 1 for classes in representations.values())
    # Explicitly include a non-palindromic RC pair and a non-equivalent pair.
    assert np.array_equal(canonical_one_hot("ACGTCAG"), canonical_one_hot(reverse_complement("ACGTCAG")))
    assert not np.array_equal(canonical_one_hot("ACGTCAG"), canonical_one_hot("ACGTTAG"))


def test_legacy_average_has_cross_class_collision():
    sequences = all_kmers()
    representations = {}
    for sequence in sequences:
        representations.setdefault(tuple(rc_symmetric_one_hot(sequence)), set()).add(canonical_rc(sequence))
    assert sum(len(classes) > 1 for classes in representations.values()) == 1872


def test_all_v06_model_scores_are_rc_invariant_for_all_oriented_7mers():
    candidates = all_kmers()
    protein = torch.zeros(4)
    target = "TAGCAGGATGTGT"
    models = [
        (CandidateDNAOnlyV06(8), CandidateDNAOnlyV06.model_name),
        (ProteinCandidateV06(4, 8), ProteinCandidateV06.model_name),
        (TargetCandidateOnlyV06(8), TargetCandidateOnlyV06.model_name),
        (ProteinTargetCandidateV06(4, 8), ProteinTargetCandidateV06.model_name),
    ]
    rc_candidates = [reverse_complement(sequence) for sequence in candidates]
    for model, name in models:
        first = score_rc_v06(model, name, candidates, protein_embedding=protein, target=target)
        second = score_rc_v06(model, name, rc_candidates, protein_embedding=protein, target=target)
        assert np.allclose(first, second, atol=1e-6, rtol=1e-6)


def test_target_windows_preserve_order_repeats_and_register():
    target = "AAAAAAA"
    windows = ordered_target_windows(target)
    assert len(windows) == 1
    target = "AAAAAAAA"
    windows = ordered_target_windows(target)
    assert len(windows) == 2
    assert [index for _window, index, _position in windows] == [0, 1]
    assert [position for _window, _index, position in windows] == [-1.0, 1.0]
    assert target_window_features_v06(target).shape == (2, 29)
    # A repeated window is retained rather than set-deduplicated.
    assert len(ordered_target_windows("AAAAAAAAA")) == 3


def test_target_orientation_is_invariant_after_whole_target_canonicalization():
    target = "TAGCAGGATGTGT"
    assert torch.equal(target_window_features_v06(target), target_window_features_v06(reverse_complement(target)))


def test_pooling_modes_are_available_and_position_aware():
    candidate = batch_canonical_one_hot(["ACGTCAG", "TTTTTTT"])
    target = target_window_features_v06("TAGCAGGATGTGT").unsqueeze(0).expand(2, -1, -1)
    for mode in ("mean", "max", "logsumexp", "softmax", "position_aware"):
        model = TargetCandidateOnlyV06(8, target_pooling=mode)
        output = model(target, candidate)
        assert output.shape == (2,)
        assert torch.isfinite(output).all()


def test_synthetic_conditionality_is_a_regression_gate():
    result = synthetic_conditionality_benchmark(epochs=120)
    assert len(result) == 4
    assert result["ranking_correct"].all()
    assert result["protein_swap_changes_ranking"].all()
    assert result["target_swap_changes_ranking"].all()
    assert result["protein_shuffle_not_invariant"].all()
    assert result["target_shuffle_not_invariant"].all()
    assert result["prediction_variance"].gt(0).all()


def test_score_rc_invariance_audit_passes():
    result = check_score_rc_invariance()
    assert result["score_rc_invariant"].all()
    assert result["max_abs_rc_score_difference"].le(1e-6).all()


def test_pre_replay_gates_pass():
    gates = assert_pre_replay_gates()
    assert gates["passed"].all()
