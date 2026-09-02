from pathlib import Path, PureWindowsPath

import pandas as pd

from src.natural_pbm import EXPECTED_8MER_RC_CLASSES
from src.sequence_equivalence import canonical_rc, reverse_complement


ROOT = Path(__file__).resolve().parents[1]


def test_each_training_protein_has_complete_rc_class_coverage():
    bench = pd.read_parquet(ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    counts = bench.groupby("protein_id")["canonical_rc"].nunique()
    assert counts.eq(EXPECTED_8MER_RC_CLASSES).all()
    assert not bench.duplicated(["protein_id", "canonical_rc"]).any()


def test_natural_rc_canonicalization_is_stable():
    bench = pd.read_parquet(ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    sample = bench["canonical_rc"].drop_duplicates().head(1000)
    assert sample.map(canonical_rc).eq(sample).all()
    assert sample.map(lambda seq: canonical_rc(reverse_complement(seq))).eq(sample).all()


def test_manifest_uses_project_relative_paths():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "natural_pbm_files.csv")
    paths = manifest["local_relative_path"].dropna().astype(str)
    assert not paths.empty
    for path in paths:
        assert not Path(path).is_absolute()
        assert not PureWindowsPath(path).is_absolute()
