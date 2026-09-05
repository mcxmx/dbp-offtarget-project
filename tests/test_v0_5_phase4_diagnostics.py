from pathlib import Path

import numpy as np
import pandas as pd

from src.v0_5_phase4_diagnostics import (
    DESIGNED_IDS,
    deterministic_protein_permutation,
    deterministic_target_permutation,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "v0_5"


def test_hard_case_denominators_are_consistent():
    summary = pd.read_csv(RESULTS / "hard_case_model_summary.csv")
    assert set(summary["model"]) == {"M0", "M1", "M1c", "M2", "M3"}
    assert summary["eligible"].eq(1515).all()
    assert summary["resolved"].add(summary["unresolved"]).eq(summary["eligible"]).all()
    assert summary["not_evaluable"].eq(0).all()


def test_unique_win_classification_is_nested_and_disjoint_from_common_failure():
    def keys(filename):
        frame = pd.read_csv(RESULTS / "hard_cases" / filename)
        return set(zip(frame["dbp_id"], frame["candidate"]))

    m1 = keys("m1_fail_m3_success.csv")
    m1c = keys("m1c_fail_m3_success.csv")
    m2 = keys("m2_fail_m3_success.csv")
    joint = keys("joint_controls_fail_m3_success.csv")
    all_fail = keys("all_current_models_fail.csv")
    assert joint <= m1c
    assert joint <= m2
    assert not (joint & all_fail)
    assert not (m1 & all_fail)


def test_shuffled_protein_permutation_is_deterministic_and_deranges():
    first = deterministic_protein_permutation(DESIGNED_IDS)
    second = deterministic_protein_permutation(DESIGNED_IDS)
    assert first == second
    assert set(first) == set(first.values()) == set(DESIGNED_IDS)
    assert all(left != right for left, right in first.items())


def test_shuffled_target_permutation_is_deterministic_and_deranges():
    first = deterministic_target_permutation(DESIGNED_IDS)
    second = deterministic_target_permutation(DESIGNED_IDS)
    assert first == second
    assert set(first) == set(first.values()) == set(DESIGNED_IDS)
    assert all(left != right for left, right in first.items())


def test_shuffle_diagnostics_use_frozen_models_without_retraining():
    for filename in ["shuffled_protein_diagnostic.csv", "shuffled_target_diagnostic.csv"]:
        diagnostic = pd.read_csv(RESULTS / filename)
        assert diagnostic["retrained"].eq(False).all()
        assert diagnostic["used_for_model_selection"].eq(False).all()
        expected_condition = "shuffled_protein" if "protein" in filename else "shuffled_target"
        assert set(diagnostic["condition"]) == {"original", expected_condition}


def test_pair_sampling_coverage_is_diagnostic_only_and_bounded():
    coverage = pd.read_csv(RESULTS / "pair_sampling_coverage.csv")
    assert len(coverage) == 63
    assert coverage["pairs_sampled"].eq(512).all()
    assert coverage["unique_candidates_in_pairs"].between(1, 8192).all()
    assert coverage["candidate_coverage_fraction"].between(0, 1).all()
    assert coverage["n_decile_bins_covered"].between(1, 10).all()


def test_train_and_test_performance_are_explicitly_separated():
    diagnostic = pd.read_csv(RESULTS / "train_vs_test_performance.csv")
    assert set(diagnostic["partition"]) == {"train", "test"}
    assert diagnostic["diagnostic_only"].eq(True).all()
    assert diagnostic["used_for_model_selection"].eq(False).all()
    assert diagnostic[["fold_id", "seed", "dbp_id", "model", "partition"]].duplicated().sum() == 0


def test_primary_result_hashes_remain_unchanged():
    manifest = RESULTS / "PRIMARY_RESULTS_FROZEN_MANIFEST.txt"
    expected = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if line.startswith("results/"):
            relative, digest = line.split()
            expected[relative] = digest
    assert expected
    for relative, digest in expected.items():
        assert sha256_file(ROOT / relative) == digest


def test_phase4_replay_matches_frozen_primary():
    replay = pd.read_csv(RESULTS / "phase4_replay_validation.csv")
    assert len(replay) == 105
    assert replay["matches_frozen_primary"].eq(True).all()
    assert np.isclose(replay["absolute_difference"].max(), 0.0, atol=1e-10)
