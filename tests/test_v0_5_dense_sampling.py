from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src.v0_5_dense_sampling import audit_pair_sampling, sample_dense_pairs
from src.v0_5_dense_training import select_untouched_smoke_fold
from src.v0_5_models import CandidateDNAOnly, ProteinCandidate, ProteinTargetCandidate, TargetCandidateOnly
from src.v0_5_training import sample_rank_pairs


ROOT = Path(__file__).resolve().parents[1]


def _toy_benchmark(n: int = 100) -> pd.DataFrame:
    bases = "ACGT"
    sequences = []
    for index in range(n):
        value = index
        chars = []
        for _ in range(7):
            chars.append(bases[value % 4])
            value //= 4
        sequences.append("".join(chars))
    return pd.DataFrame(
        {
            "protein_id": ["P1"] * n,
            "candidate_dna": sequences,
            "experimental_score": np.linspace(-1.0, 1.0, n),
        }
    )


def test_dense_sampler_is_deterministic_and_exact():
    benchmark = _toy_benchmark()
    first = sample_dense_pairs(
        benchmark,
        ["P1"],
        pairs_per_protein=200,
        seed=17,
        protocol="D4096",
    )
    second = sample_dense_pairs(
        benchmark,
        ["P1"],
        pairs_per_protein=200,
        seed=17,
        protocol="D4096",
    )
    pd.testing.assert_frame_equal(first, second)
    assert len(first) == 200
    assert first["difficulty"].value_counts().to_dict() == {"easy": 80, "medium": 70, "hard": 50}


def test_dense_sampler_is_within_protein_and_has_no_unordered_duplicates():
    pairs = sample_dense_pairs(_toy_benchmark(), ["P1"], pairs_per_protein=200, seed=17, protocol="D4096")
    keys = {tuple(sorted((int(row.left_index), int(row.right_index)))) for row in pairs.itertuples()}
    assert len(keys) == len(pairs)
    assert set(pairs["protein_id"]) == {"P1"}


def test_real_sampling_audit_meets_registered_coverage_targets():
    audit = pd.read_csv(ROOT / "results" / "v0_5_dense" / "pair_sampling_audit.csv")
    for protocol, threshold in [("S512", 0.0), ("D4096", 0.60), ("D16384", 0.90)]:
        subset = audit.loc[audit["protocol"].eq(protocol)]
        assert len(subset) == 7
        assert (subset["candidate_coverage"] >= threshold).all()
        assert subset["n_rank_deciles_covered"].eq(10).all()
        assert subset["duplicate_pair_fraction"].eq(0.0).all()
    assert audit.loc[audit["protocol"].eq("S512"), "pair_count"].eq(512).all()
    assert audit.loc[audit["protocol"].eq("D4096"), "pair_count"].eq(4096).all()
    assert audit.loc[audit["protocol"].eq("D16384"), "pair_count"].eq(16384).all()


def test_dense_sampling_does_not_change_benchmark_labels():
    benchmark = _toy_benchmark()
    before = benchmark.copy(deep=True)
    sample_dense_pairs(benchmark, ["P1"], pairs_per_protein=200, seed=17, protocol="D4096")
    pd.testing.assert_frame_equal(benchmark, before)


def test_s512_sampler_remains_the_frozen_sampler():
    benchmark = _toy_benchmark()
    first = sample_rank_pairs(benchmark, ["P1"], pairs_per_protein=20, seed=42)
    second = sample_rank_pairs(benchmark, ["P1"], pairs_per_protein=20, seed=42)
    pd.testing.assert_frame_equal(first, second)


def test_dense_config_models_are_the_frozen_global_family():
    config = __import__("json").loads(
        (ROOT / "metadata" / "v0_5_dense" / "dense_training_config.json").read_text(encoding="utf-8")
    )
    assert config["models"] == [
        CandidateDNAOnly.model_name,
        "M1c_ProteinCandidateCapacityMatched",
        TargetCandidateOnly.model_name,
        ProteinTargetCandidate.model_name,
    ]
    assert config["epochs"] == 18
    assert config["seed"] == 17


def test_first_smoke_fold_is_entirely_still_untouched():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_5" / "v0_5_split_manifest.csv")
    fold_id, train, test = select_untouched_smoke_fold(
        manifest,
        still_untouched={"DBP48", "DBP6", "DBP9"},
    )
    assert fold_id == "protein_cluster_loco_fold_3"
    assert test == ["DBP48"]
    assert "DBP48" not in train


def test_dense_protocol_does_not_change_model_architecture_contract():
    assert CandidateDNAOnly.model_name == "M0_CandidateDNAOnly"
    assert ProteinCandidate.model_name == "M1_ProteinCandidate"
    assert TargetCandidateOnly.model_name == "M2_TargetCandidateOnly"
    assert ProteinTargetCandidate.model_name == "M3_ProteinTargetCandidate"
