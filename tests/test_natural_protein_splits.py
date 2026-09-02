from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_cluster_40_does_not_cross_splits():
    splits = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "natural_pbm_splits.csv")
    assert {"train", "validation", "natural_test"}.issubset(set(splits["split"]))
    assert splits.groupby("cluster_40")["split"].nunique().le(1).all()


def test_split_counts_are_nontrivial():
    splits = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "natural_pbm_splits.csv")
    counts = splits["split"].value_counts()
    assert counts["train"] >= 30
    assert counts["validation"] >= 5
    assert counts["natural_test"] >= 5


def test_training_proteins_have_sequences_and_no_designed_ids():
    bench = pd.read_parquet(ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    splits = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "natural_pbm_splits.csv")
    train_ids = set(splits.loc[splits["split"] == "train", "protein_id"])
    train = bench[bench["protein_id"].isin(train_ids)]
    assert train["protein_sequence"].notna().all()
    assert not train["protein_id"].str.contains("DBP", regex=False).any()
