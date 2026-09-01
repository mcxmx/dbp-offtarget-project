from pathlib import Path

import pandas as pd
from pandas.api.types import is_numeric_dtype

from src.sequence_equivalence import canonical_rc


ROOT = Path(__file__).resolve().parents[1]
ALLOWED_PROTEINS = {"DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"}
EXPECTED_ORIENTED = 4**7
EXPECTED_RC_CLASSES = EXPECTED_ORIENTED // 2


def test_upbm_long_dna_and_score_qc():
    long_df = pd.read_parquet(ROOT / "data" / "interim" / "gse237017" / "upbm_7mers_long.parquet")
    assert set(long_df["protein_id"]) == ALLOWED_PROTEINS
    assert long_df["gsm_id"].notna().all()
    assert long_df["dna_7mer"].str.len().eq(7).all()
    assert long_df["dna_7mer"].str.fullmatch("[ACGT]{7}").all()
    for column in ["e_score", "median_intensity", "z_score"]:
        assert column in long_df.columns
        assert is_numeric_dtype(long_df[column])


def test_each_sample_has_oriented_sequences_and_rc_classes():
    long_df = pd.read_parquet(ROOT / "data" / "interim" / "gse237017" / "upbm_7mers_long.parquet")
    long_df = long_df.copy()
    long_df["canonical_7mer"] = long_df["dna_7mer"].map(canonical_rc)
    per_sample = long_df.groupby("gsm_id").agg(
        n_rows=("dna_7mer", "size"),
        n_oriented=("dna_7mer", "nunique"),
        n_rc_classes=("canonical_7mer", "nunique"),
    )
    assert per_sample["n_rows"].eq(EXPECTED_ORIENTED).all()
    assert per_sample["n_oriented"].eq(EXPECTED_ORIENTED).all()
    assert per_sample["n_rc_classes"].eq(EXPECTED_RC_CLASSES).all()


def test_paired_reverse_complement_scores_are_consistent():
    long_df = pd.read_parquet(ROOT / "data" / "interim" / "gse237017" / "upbm_7mers_long.parquet")
    long_df = long_df.copy()
    long_df["canonical_7mer"] = long_df["dna_7mer"].map(canonical_rc)
    diff = (
        long_df.groupby(["gsm_id", "canonical_7mer"])["e_score"]
        .agg(lambda values: float(values.max() - values.min()))
        .max()
    )
    assert diff == 0.0


def test_manifest_hash_and_sample_metadata_present():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_3" / "gse237017_file_manifest.csv")
    downloaded = manifest[manifest["download_status"] == "downloaded"]
    assert not downloaded.empty
    assert downloaded["sha256"].str.fullmatch("[0-9a-f]{64}").all()
    samples = pd.read_csv(ROOT / "metadata" / "v0_3" / "gse237017_samples.csv", dtype={"replicate": str})
    assert samples["gsm_id"].notna().all()
    assert samples["source_url"].str.startswith("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM").all()
    assert not samples.duplicated(["protein_id", "replicate"]).any()


def test_sample_coverage_qc_matches_expected_counts():
    coverage = pd.read_csv(ROOT / "results" / "v0_3" / "tables" / "sample_coverage_qc.csv")
    assert coverage["n_unique_7mers"].eq(EXPECTED_ORIENTED).all()
    assert coverage["n_missing_7mers"].eq(0).all()
    assert coverage["n_duplicate_7mers"].eq(0).all()
    assert coverage["fraction_complete"].eq(1.0).all()
