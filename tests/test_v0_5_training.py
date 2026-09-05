from pathlib import Path

import numpy as np
import pandas as pd
import torch

from src.v0_5_training import V05Config, pairwise_ranking_loss, sample_rank_pairs
from src.v0_5_models import CandidateDNAOnly, ProteinCandidate, ProteinTargetCandidate, model_parameter_counts


ROOT = Path(__file__).resolve().parents[1]


def _toy_benchmark() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein_id": ["P1"] * 8 + ["P2"] * 8,
            "candidate_dna": [
                "AAAAAAA",
                "AAAAAAC",
                "AAAAAAG",
                "AAAAAAT",
                "CCCCCCC",
                "GGGGGGG",
                "TTTTTTT",
                "ACGTACG",
            ]
            * 2,
            "experimental_score": list(np.linspace(0.0, 1.0, 8)) * 2,
        }
    )


def test_pairwise_ranking_loss_has_correct_sign():
    good = pairwise_ranking_loss(
        torch.tensor([2.0]),
        torch.tensor([0.0]),
        torch.tensor([1.0]),
        torch.tensor([0.0]),
    )
    bad = pairwise_ranking_loss(
        torch.tensor([0.0]),
        torch.tensor([2.0]),
        torch.tensor([1.0]),
        torch.tensor([0.0]),
    )
    assert good.item() < bad.item()


def test_pair_sampler_is_deterministic_and_within_protein():
    benchmark = _toy_benchmark()
    first = sample_rank_pairs(benchmark, ["P1"], pairs_per_protein=20, seed=42)
    second = sample_rank_pairs(benchmark, ["P1"], pairs_per_protein=20, seed=42)
    pd.testing.assert_frame_equal(first, second)
    assert set(first["protein_id"]) == {"P1"}


def test_pair_sampler_excludes_held_out_proteins_by_input_contract():
    benchmark = _toy_benchmark()
    train_proteins = ["P1"]
    pairs = sample_rank_pairs(benchmark, train_proteins, pairs_per_protein=20, seed=42)
    assert not set(pairs["protein_id"]) - set(train_proteins)


def test_smoke_fold_is_defined_by_split_manifest():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_5" / "v0_5_split_manifest.csv")
    fold = manifest.loc[
        manifest["split_name"].eq("protein_cluster_loco")
        & manifest["fold_id"].eq("protein_cluster_loco_fold_1")
    ]
    assert set(fold.loc[fold["partition"].eq("test"), "dbp_id"]) == {"DBP1", "DBP3"}
    assert not (
        set(fold.loc[fold["partition"].eq("train"), "dbp_id"])
        & set(fold.loc[fold["partition"].eq("test"), "dbp_id"])
    )


def test_phase2_config_is_fixed_and_not_test_tuned():
    config = V05Config()
    assert config.seed == 42
    assert config.smoke_fold_id == "protein_cluster_loco_fold_1"
    assert config.epochs == 18
    assert config.pair_count_per_protein == 512


def test_m1_and_m3_have_comparable_parameter_counts():
    m1 = model_parameter_counts(ProteinCandidate(480, hidden_dim=32))[0]
    m3 = model_parameter_counts(ProteinTargetCandidate(480, hidden_dim=32))[0]
    assert m3 / m1 < 2.0
    assert m3 / m1 > 0.5


def test_m0_has_no_protein_parameters():
    trainable, frozen = model_parameter_counts(CandidateDNAOnly(hidden_dim=32))
    assert trainable > 0
    assert frozen == 0
