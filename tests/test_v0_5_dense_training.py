from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DENSE_RESULTS = ROOT / "results" / "v0_5_dense"
PRIMARY_RESULTS = ROOT / "results" / "v0_5"


def test_dense_smoke_uses_only_training_proteins_for_pairs():
    info = pd.read_csv(DENSE_RESULTS / "smoke_split_info.csv").iloc[0]
    training = set(str(info["train_proteins"]).split("|"))
    held_out = set(str(info["test_proteins"]).split("|"))
    pairs = pd.read_csv(DENSE_RESULTS / "smoke_training_pairs.csv")
    assert training.isdisjoint(held_out)
    assert set(pairs["protein_id"]).issubset(training)
    assert not set(pairs["protein_id"]).intersection(held_out)


def test_dense_smoke_has_exact_protocol_pair_counts_and_rc_units():
    pairs = pd.read_csv(DENSE_RESULTS / "smoke_training_pairs.csv")
    health = pd.read_csv(DENSE_RESULTS / "training_health.csv")
    expected_per_protein = {"S512": 512, "D4096": 4096, "D16384": 16384}
    for protocol, expected in expected_per_protein.items():
        counts = pairs.loc[pairs["protocol"].eq(protocol)].groupby("protein_id").size()
        assert counts.nunique() == 1
        assert counts.iloc[0] == expected
    evaluation = pd.read_csv(DENSE_RESULTS / "smoke_results.csv")
    assert evaluation["n_rc_units"].eq(8192).all()
    assert health["training_pair_count"].isin({3072, 24576, 98304}).all()


def test_dense_smoke_s512_replays_frozen_primary_fold3_seed17():
    dense = pd.read_csv(DENSE_RESULTS / "smoke_results.csv")
    dense = dense.loc[
        dense["protocol"].eq("S512")
        & dense["partition"].eq("test")
        & dense["dbp_id"].eq("DBP48")
    ].set_index("model")["spearman"]
    primary = pd.read_csv(PRIMARY_RESULTS / "primary_seed_level_results.csv")
    primary = primary.loc[
        primary["fold_id"].eq("protein_cluster_loco_fold_3")
        & primary["seed"].eq(17)
        & primary["dbp_id"].eq("DBP48")
    ].set_index("model")["spearman"]
    for model in ["M0", "M1c", "M2", "M3"]:
        assert np.isclose(dense[model], primary[model], atol=1e-10, rtol=1e-10)


def test_dense_smoke_shuffle_diagnostics_are_inference_only():
    shuffle = pd.read_csv(DENSE_RESULTS / "shuffle_diagnostics.csv")
    assert shuffle["retrained"].eq(False).all()
    assert set(shuffle["protocol"]) == {"S512", "D4096", "D16384"}
    assert set(shuffle["model"]) == {"M0", "M1c", "M2", "M3"}
    assert shuffle.loc[shuffle["model"].ne("M0"), "prediction_correlation"].notna().all()


def test_dense_smoke_exposure_bookkeeping_is_frozen():
    info = pd.read_csv(DENSE_RESULTS / "smoke_split_info.csv").iloc[0]
    manifest = pd.read_csv(ROOT / "metadata" / "v0_5_dense" / "development_exposure_manifest.csv")
    assert info["fold_id"] == "protein_cluster_loco_fold_3"
    assert info["test_proteins"] == "DBP48"
    status = manifest.set_index("protein_id")["phase6a_exposure_status"].to_dict()
    assert status["DBP48"] == "development_exposed_for_dense_supervision"
    assert status["DBP6"] == "still_untouched"
    assert status["DBP9"] == "still_untouched"


def test_dense_smoke_health_has_no_technical_failures():
    health = pd.read_csv(DENSE_RESULTS / "training_health.csv")
    assert health["status"].eq("complete").all()
    assert health["nan_inf_count"].eq(0).all()
    assert health["test_prediction_variance"].gt(0).all()
    assert health[["first_epoch_loss", "final_epoch_loss"]].notna().all().all()
