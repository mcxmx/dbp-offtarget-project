from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_deeppbs_provenance_records_official_repo_and_commit():
    text = (ROOT / "external" / "deeppbs" / "PROVENANCE.md").read_text(encoding="utf-8")
    assert "https://github.com/timkartar/DeepPBS" in text
    assert "8bfb211dd67f02877841f6f33aa493ddf7daedf9" in text
    assert "Docker/WSL" in text


def test_deeppbs_prediction_schema_empty_when_not_run():
    pred = pd.read_parquet(ROOT / "results" / "v0_4_1" / "tables" / "deeppbs_designed_predictions.parquet")
    required = {
        "protein_id",
        "canonical_7mer",
        "deeppbs_score",
        "structure_id",
        "structure_type",
        "model_version",
        "overlap_status",
    }
    assert required.issubset(pred.columns)
    assert pred.empty


def test_deeppbs_overlap_and_weight_manifests_exist():
    overlap = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "deeppbs_overlap_audit.csv")
    weights = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "deeppbs_weight_manifest.csv")
    assert set(overlap["protein_id"]) == {"DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"}
    assert weights["sha256"].str.fullmatch("[0-9a-f]{64}").all()
