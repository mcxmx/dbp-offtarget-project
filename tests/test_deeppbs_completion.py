from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.deeppbs_adapter import (
    parse_prediction_npz,
    pwm_oriented_score,
    pwm_rc_class_score,
)
from src.sequence_equivalence import canonical_rc, reverse_complement


ROOT = Path(__file__).resolve().parents[1]


def test_official_prediction_npz_parses_with_fixed_base_order():
    path = (
        ROOT
        / "results"
        / "v0_4_2"
        / "external_runs"
        / "deeppbs_official_example"
        / "5x6g.npz_predict.npz"
    )
    parsed = parse_prediction_npz(path)
    assert parsed["base_order"] == ("A", "C", "G", "T")
    assert parsed["length"] == 14
    assert parsed["sequence"] == "TCAGTCTAGACATA"


def test_prediction_parser_rejects_wrong_probability_rows(tmp_path):
    path = tmp_path / "bad.npz"
    np.savez(path, P=np.ones((7, 4)), Seq=np.eye(4)[[0, 1, 2, 3, 0, 1, 2]])
    with pytest.raises(ValueError, match="not normalized"):
        parse_prediction_npz(path)


def test_pwm_scorer_has_known_log_probability_value():
    pwm = np.full((7, 4), 0.0)
    pwm[:, 0] = 0.9
    pwm[:, 1:] = 0.1 / 3.0
    score, offset = pwm_oriented_score(pwm, "AAAAAAA")
    assert offset == 0
    assert np.isclose(score, 7.0 * np.log(0.9))


def test_pwm_rc_class_score_is_orientation_invariant():
    rng = np.random.default_rng(42)
    pwm = rng.random((9, 4))
    pwm /= pwm.sum(axis=1, keepdims=True)
    sequence = "ACGTTGC"
    canonical = canonical_rc(sequence)
    result = pwm_rc_class_score(pwm, canonical)
    assert result["canonical_7mer"] == canonical
    assert np.isclose(
        result["prediction_score"],
        max(
            pwm_oriented_score(pwm, canonical)[0],
            pwm_oriented_score(pwm, reverse_complement(canonical))[0],
        ),
    )
    assert canonical_rc(reverse_complement(canonical)) == canonical


def test_pwm_scorer_rejects_profiles_shorter_than_seven():
    with pytest.raises(ValueError, match="at least 7"):
        pwm_oriented_score(np.full((6, 4), 0.25), "AAAAAAA")


def test_completed_deeppbs_landscape_has_unique_rc_class_units():
    path = (
        ROOT
        / "results"
        / "v0_4_2"
        / "tables"
        / "deeppbs_predictions_completed_v0_4_2.parquet"
    )
    table = pd.read_parquet(path)
    assert len(table) == 2 * 8192
    assert set(table["protein_id"]) == {"DBP35", "DBP48"}
    assert table.duplicated(["protein_id", "canonical_7mer"]).sum() == 0
    assert table.groupby("protein_id").size().to_dict() == {
        "DBP35": 8192,
        "DBP48": 8192,
    }
    assert table["deeppbs_score"].notna().all()
    assert table["experimental_E_score"].notna().all()
    assert (
        table["experimental_score_type"]
        .eq("processed experimental uPBM E-score consensus")
        .all()
    )


def test_completed_performance_keeps_missing_predictions_missing():
    path = (
        ROOT
        / "results"
        / "v0_4_2"
        / "tables"
        / "deeppbs_performance_completed_v0_4_2.csv"
    )
    performance = pd.read_csv(path).set_index("protein_id")
    assert performance.loc["DBP35", "status"] == "evaluated"
    assert performance.loc["DBP48", "status"] == "evaluated"
    assert performance.loc["DBP35", "n_rc_classes"] == 8192
    assert performance.loc["DBP48", "n_rc_classes"] == 8192
    for protein_id in ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9"]:
        assert performance.loc[protein_id, "status"] == "not_evaluable_missing_prediction"
        assert pd.isna(performance.loc[protein_id, "spearman"])


def test_completed_run_manifest_has_correct_contacts_and_portable_paths():
    manifest = pd.read_csv(
        ROOT
        / "results"
        / "v0_4_2"
        / "tables"
        / "deeppbs_run_manifest_completed_v0_4_2.csv"
    ).set_index("protein_id")
    assert manifest.loc["DBP35", "contact_count"] == 266
    assert manifest.loc["DBP48", "contact_count"] == 224
    text = manifest.to_csv()
    assert "E:\\" not in text
    assert "C:\\" not in text
    assert (
        ROOT
        / "results"
        / "v0_4_2"
        / "external_runs"
        / "deeppbs_official_example"
        / "official_example_status.txt"
    ).read_text(encoding="utf-8").startswith("OFFICIAL_EXAMPLE_PASS")
