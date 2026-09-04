from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DESIGNED = {"DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"}


def test_disagreement_resolution_preserves_candidate_denominator():
    resolution = pd.read_csv(
        ROOT
        / "results"
        / "v0_4_2"
        / "tables"
        / "disagreement_resolution_v0_4_2.csv"
    )
    assert set(resolution["protein_id"]) == DESIGNED
    assert set(resolution["method"]) == {
        "sequence_kmer3",
        "SimpleProteinConditionalBaseline",
        "FrozenPLMProteinConditionalBaseline",
        "NA-MPNN diagnostic",
        "DeepPBS",
    }
    assert resolution.groupby("method")["n_total_candidates"].sum().eq(1515).all()
    assert resolution.loc[
        resolution["method"].eq("DeepPBS"), "evaluation_status"
    ].eq("not_evaluable_missing_prediction").all()
    assert resolution.loc[
        resolution["method"].eq("DeepPBS"), "n_evaluable"
    ].eq(0).all()


def test_difficulty_factors_have_training_space_and_performance_columns():
    factors = pd.read_csv(
        ROOT / "results" / "v0_4_2" / "tables" / "designed_difficulty_factors.csv"
    )
    assert set(factors["protein_id"]) == DESIGNED
    assert factors[
        [
            "max_natural_train_sequence_identity",
            "nearest_natural_train_esm_euclidean_distance",
            "sequence_kmer3_spearman",
            "simple_pc_spearman",
            "frozen_plm_spearman",
        ]
    ].notna().all().all()


def test_v0_4_2_diagnostic_summary_reports_rc_class_rows():
    summary = pd.read_json(
        ROOT / "results" / "v0_4_2" / "tables" / "v0_4_2_diagnostic_summary.json"
    )
    assert int(summary.loc[0, "n_designed_rows"]) == 57344
    assert int(summary.loc[0, "n_disagreement_candidates"]) == 1515
    assert int(summary.loc[0, "n_common_high_experiment_low_all_core"]) == 263
