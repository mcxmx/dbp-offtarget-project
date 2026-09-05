from pathlib import Path

import numpy as np
import pandas as pd

from src.v0_5_primary_evaluation import (
    MODEL_ORDER,
    UNSEEN_PRIMARY,
    aggregate_seed_level,
    build_per_protein_results,
    build_primary_macro_summary,
    validate_loco_manifest,
)
from src.v0_5_training import V05Config, run_smoke


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "metadata" / "v0_5"


def test_all_primary_folds_cover_each_dbp_once_without_cluster_leakage():
    splits = pd.read_csv(METADATA / "v0_5_split_manifest.csv")
    audit = validate_loco_manifest(splits, "protein_cluster_loco", "protein_cluster")
    assert audit["n_folds"] == 4
    assert sorted(audit["all_proteins"]) == ["DBP1", "DBP3", "DBP35", "DBP48", "DBP5", "DBP6", "DBP9"]


def test_all_strict_folds_cover_each_dbp_once_without_component_leakage():
    splits = pd.read_csv(METADATA / "v0_5_split_manifest.csv")
    audit = validate_loco_manifest(splits, "combined_component_loco", "combined_component")
    assert audit["n_folds"] == 3
    assert len(audit["all_proteins"]) == 7


def test_seed_aggregation_and_macro_median_are_per_protein_first():
    rows = []
    for dbp, value in [("DBP1", 0.1), ("DBP3", 0.3)]:
        for seed, offset in [(17, 0.0), (29, 0.1), (43, -0.1)]:
            rows.append(
                {
                    "split_type": "protein_cluster_loco",
                    "fold_id": "fold1",
                    "seed": seed,
                    "dbp_id": dbp,
                    "protein_cluster": "cluster1",
                    "combined_component": "component1",
                    "model": "M0",
                    "n_units": 8192,
                    "spearman": value + offset,
                    "ndcg_1pct": 0.5,
                    "ndcg_5pct": 0.5,
                    "pairwise_accuracy": 0.5,
                    "top1pct_recovery": 0.1,
                    "training_proteins": "DBP5",
                    "test_proteins": "DBP1|DBP3",
                    "runtime_seconds": 1.0,
                    "training_pair_count": 512,
                    "development_exposed": dbp in {"DBP1", "DBP3"},
                    "status": "complete",
                }
            )
    seed_level = pd.DataFrame(rows)
    aggregated = aggregate_seed_level(seed_level)
    assert len(aggregated) == 2
    assert np.isclose(aggregated["spearman"].mean(), 0.2)
    per_protein = build_per_protein_results(seed_level)
    assert np.allclose(per_protein["M0"], [0.1, 0.3])
    macro = build_primary_macro_summary(per_protein, (17, 29, 43))
    assert set(macro["model"]) == {"M0", *MODEL_ORDER[1:]}
    assert np.isclose(macro.loc[macro["model"].eq("M0"), "all7_macro_median"].iloc[0], 0.2)


def test_development_exposure_and_unseen_sets_are_frozen():
    assert UNSEEN_PRIMARY == {"DBP5", "DBP35", "DBP48", "DBP6", "DBP9"}
    config = V05Config()
    assert config.evaluation_seeds == (17, 29, 43)
    assert config.seed == 42


def test_primary_manifest_uses_intended_target_not_pbm_motif():
    manifest = pd.read_csv(METADATA / "designed_target_manifest_v0_5.csv")
    assert (manifest["primary_target"] == manifest["intended_design_target"]).all()
    assert (manifest["primary_target"] != manifest["pbm_motif"]).all()


def test_config_does_not_permit_designed_test_selection():
    config = V05Config()
    assert config.as_dict()["designed_test_used_for_selection"] is False


def test_existing_primary_artifacts_have_complete_unique_coverage():
    results = ROOT / "results" / "v0_5"
    seed_level = pd.read_csv(results / "primary_seed_level_results.csv")
    strict_seed_level = pd.read_csv(results / "strict_component_seed_level_results.csv")
    expected_models = {"M0", "M1", "M1c", "M2", "M3"}
    assert set(seed_level["model"]) == expected_models
    assert set(strict_seed_level["model"]) == expected_models
    assert seed_level[["fold_id", "seed", "dbp_id", "model"]].duplicated().sum() == 0
    assert strict_seed_level[["fold_id", "seed", "dbp_id", "model"]].duplicated().sum() == 0
    assert seed_level["n_units"].eq(8192).all()
    assert strict_seed_level["n_units"].eq(8192).all()
    assert seed_level["status"].eq("complete").all()
    assert strict_seed_level["status"].eq("complete").all()
    assert len(seed_level) == 105
    assert len(strict_seed_level) == 105


def test_primary_training_health_has_no_technical_failures():
    health = pd.read_csv(ROOT / "results" / "v0_5" / "primary_training_health.csv")
    assert health["nan_inf_count"].eq(0).all()
    assert health["prediction_variance"].gt(0).all()
    assert health["status"].eq("complete").all()


def test_exposure_flags_match_the_development_note():
    per_protein = pd.read_csv(ROOT / "results" / "v0_5" / "primary_per_protein_results.csv")
    exposed = set(per_protein.loc[per_protein["development_exposed"], "dbp_id"])
    assert exposed == {"DBP1", "DBP3"}
    assert set(per_protein.loc[~per_protein["development_exposed"], "dbp_id"]) == UNSEEN_PRIMARY


def test_smoke_output_is_reproducible_under_frozen_config():
    expected = pd.read_csv(ROOT / "results" / "v0_5" / "smoke_test_results.csv")
    actual = run_smoke(V05Config())["evaluation"]
    key = ["protein_id", "model"]
    expected = expected.sort_values(key).reset_index(drop=True)
    actual = actual.sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(expected, actual, check_exact=False, atol=1e-7, rtol=1e-7)
