from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_no_exact_designed_sequence_in_natural_training_split():
    designed = pd.read_csv(ROOT / "metadata" / "v0_3" / "designed_dbp_sequences.csv")
    designed_sequences = set(designed["protein_sequence"].dropna())
    natural = pd.read_parquet(ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    splits = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "natural_pbm_splits.csv")
    train_ids = set(splits.loc[splits["split"] == "train", "protein_id"])
    natural_train_sequences = set(natural.loc[natural["protein_id"].isin(train_ids), "protein_sequence"].dropna())
    assert designed_sequences.isdisjoint(natural_train_sequences)


def test_no_absolute_windows_paths_in_v0_4_1_metadata():
    for path in (ROOT / "metadata" / "v0_4_1").glob("*.csv"):
        text = path.read_text(encoding="utf-8")
        assert "E:\\" not in text
        assert "C:\\" not in text
