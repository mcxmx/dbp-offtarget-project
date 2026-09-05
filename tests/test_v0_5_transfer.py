from hashlib import sha256
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_DESIGNED = {"DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"}


def test_current_exposure_manifest_marks_all_designed_proteins_exposed():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_5_transfer" / "exposure_manifest.csv")
    assert set(manifest["protein_id"]) == EXPECTED_DESIGNED
    assert manifest["current_exposure_status"].eq("development_exposed").all()
    assert manifest["confirmatory_status"].eq("not_available").all()


def test_transfer_training_regime_manifest_records_stop_and_no_designed_labels_in_natural():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_5_transfer" / "training_regime_manifest.csv")
    natural = manifest.loc[manifest["regime"].eq("R-NATURAL")].iloc[0]
    assert natural["designed_labels_used"] == "no"
    assert natural["bridge_status"] == "not_run_after_explicit_h3_stop"
    designed = manifest.loc[manifest["regime"].eq("R-DESIGNED")].iloc[0]
    assert designed["designed_fine_tuning_allowed"] == "fold-training proteins only"
    transfer = manifest.loc[manifest["regime"].eq("R-NATURAL+DESIGNED")].iloc[0]
    assert transfer["held_out_designed_cluster_excluded"] == "yes"
    assert transfer["target_features_used"] == "no"


def test_transfer_config_records_no_target_features_and_no_designed_labels_in_natural_training():
    config = pd.read_json(ROOT / "metadata" / "v0_5_transfer" / "transfer_config.json", typ="series")
    assert bool(config["designed_labels_added_to_natural_training"]) is False
    assert bool(config["target_conditioned_model_implemented"]) is False
    assert bool(config["target_features_used_in_natural_training"]) is False


def test_historical_dense_exposure_snapshot_is_not_rewritten():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_5_dense" / "development_exposure_manifest.csv")
    manifest = manifest.set_index("protein_id")
    assert manifest.loc["DBP48", "phase6a_exposure_status"] == "development_exposed_for_dense_supervision"
    assert manifest.loc["DBP6", "phase6a_exposure_status"] == "still_untouched"
    assert manifest.loc["DBP9", "phase6a_exposure_status"] == "still_untouched"


def test_frozen_v0_5_primary_hashes_are_unchanged():
    freeze = ROOT / "results" / "v0_5" / "PRIMARY_RESULTS_FROZEN_MANIFEST.txt"
    for line in freeze.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) != 2 or not parts[0].startswith("results/"):
            continue
        relative_path, expected_hash = parts
        actual_hash = sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash, relative_path


def test_natural_and_designed_benchmarks_remain_separate():
    natural = pd.read_parquet(ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    designed = pd.read_parquet(ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet")
    assert natural["dna_length"].eq(8).all()
    assert designed["canonical_7mer"].str.len().eq(7).all()
    assert set(natural["protein_id"]).isdisjoint(EXPECTED_DESIGNED)


def test_phase7a_stop_artifacts_only_record_prior_context_and_no_new_bridge():
    summary = pd.read_csv(ROOT / "results" / "v0_5_transfer" / "transfer_stop_summary.csv")
    assert summary.iloc[0]["decision"] == "NOT_SUPPORTED"
    assert summary.iloc[0]["direct_bridge_status"] == "not_run_after_explicit_stop"
    replay = pd.read_csv(ROOT / "results" / "v0_5_transfer" / "simplepc_replay.csv")
    assert replay.iloc[0]["new_replay_performed"] == False
    assert replay.iloc[0]["reproduction_quality"] == "not_assessed"


def test_phase7a_report_records_h3_stop_and_future_validation_requirement():
    report = (ROOT / "docs" / "v0_5_transfer" / "PHASE7A_TRAINING_DIVERSITY_REPORT.md").read_text(encoding="utf-8")
    assert "H3: NOT SUPPORTED" in report
    assert "independent designed-DBP dataset" in report
    assert "H4:" in report
