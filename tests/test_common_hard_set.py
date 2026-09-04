from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_common_hard_set_uses_core_rankings_without_zero_filling():
    path = ROOT / "data" / "processed" / "v0_4_2" / "common_hard_specificity_cases.parquet"
    hard = pd.read_parquet(path)
    required = {
        "protein_id",
        "canonical_7mer",
        "experimental_score",
        "experimental_percentile",
        "sequence_score_percentile",
        "simple_pc_score_percentile",
        "frozen_plm_score_percentile",
        "failure_type",
    }
    assert required.issubset(hard.columns)
    assert not hard.duplicated(["protein_id", "canonical_7mer"]).any()
    assert hard["experimental_percentile"].ge(0.95).all()
    assert hard[
        [
            "sequence_score_percentile",
            "simple_pc_score_percentile",
            "frozen_plm_score_percentile",
        ]
    ].le(0.50).all().all()
    assert hard[
        ["sequence_score", "simple_pc_score", "frozen_plm_score"]
    ].notna().all().all()
    assert hard["failure_type"].eq("high_experiment_low_all_core").all()


def test_common_hard_set_covers_all_designed_proteins():
    hard = pd.read_parquet(
        ROOT / "data" / "processed" / "v0_4_2" / "common_hard_specificity_cases.parquet"
    )
    assert set(hard["protein_id"]) == {
        "DBP1",
        "DBP3",
        "DBP5",
        "DBP6",
        "DBP9",
        "DBP35",
        "DBP48",
    }
