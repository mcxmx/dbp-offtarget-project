from pathlib import Path

import pandas as pd
from pandas.api.types import is_numeric_dtype

from src.natural_pbm import read_uniprobe_8mer_file


ROOT = Path(__file__).resolve().parents[1]


def test_uniprobe_8mer_parser_validates_core_columns():
    raw = b"ACGTACGT\tACGTACGT\t0.25\t100.0\t1.5\nAAAAAAAC\tGTTTTTTT\t-0.1\t50.0\t-0.2\n"
    df = read_uniprobe_8mer_file(raw)
    assert set(["dna_sequence", "experimental_score", "median_intensity", "z_score", "canonical_rc"]).issubset(df.columns)
    assert df["dna_sequence"].str.fullmatch("[ACGT]{8}").all()
    assert is_numeric_dtype(df["experimental_score"])


def test_natural_benchmark_schema_and_provenance():
    bench = pd.read_parquet(ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    required = {
        "protein_id",
        "protein_sequence",
        "protein_family",
        "species",
        "dna_sequence",
        "dna_length",
        "canonical_rc",
        "experimental_score",
        "experimental_percentile",
        "experiment_id",
        "assay_type",
        "score_type",
        "quality_level",
        "source",
    }
    assert required.issubset(bench.columns)
    assert bench["protein_id"].nunique() >= 50
    assert bench["dna_length"].eq(8).all()
    assert bench["dna_sequence"].str.fullmatch("[ACGT]{8}").all()
    assert is_numeric_dtype(bench["experimental_score"])
    assert bench["experiment_id"].notna().all()
    assert bench["source"].notna().all()
