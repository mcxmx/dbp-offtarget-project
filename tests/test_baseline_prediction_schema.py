from pathlib import Path

import pandas as pd
import pytest
from pandas.api.types import is_numeric_dtype

from src.models.simple_protein_conditional_baseline import SimpleProteinConditionalBaseline


ROOT = Path(__file__).resolve().parents[1]


def test_nampnn_prediction_schema_and_numeric_scores():
    pred = pd.read_parquet(ROOT / "results" / "v0_4" / "tables" / "nampnn_predictions.parquet")
    required = {
        "protein_id",
        "canonical_7mer",
        "prediction_score",
        "prediction_type",
        "structure_id",
        "model_version",
        "source_npz",
    }
    assert required.issubset(pred.columns)
    assert len(pred) == 2 * 8192
    assert pred["canonical_7mer"].str.fullmatch("[ACGT]{7}").all()
    assert is_numeric_dtype(pred["prediction_score"])
    assert pred["prediction_score"].notna().all()
    assert pred["model_version"].str.contains("NA-MPNN", regex=False).all()


def test_empty_external_prediction_files_keep_required_schema():
    for filename in ["deeppbs_predictions.parquet", "simple_pc_predictions.parquet"]:
        df = pd.read_parquet(ROOT / "data" / "processed" / "v0_4" / filename)
        required = {"protein_id", "canonical_7mer", "prediction_score", "prediction_type", "structure_id", "model_version"}
        assert required.issubset(df.columns)
        assert df.empty


def test_simple_protein_conditional_baseline_is_protein_conditioned_but_untrained():
    model = SimpleProteinConditionalBaseline()
    assert model.is_protein_conditioned is True
    assert model.training_status == "untrained_no_assay_matched_natural_training_data"
    with pytest.raises(RuntimeError):
        model.score("ACDEFGHIKLMNPQRSTVWY", "ACGTACG")
