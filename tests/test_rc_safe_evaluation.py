from pathlib import Path, PureWindowsPath

import pandas as pd

from src.sequence_equivalence import canonical_rc, has_rc_split_leakage, reverse_complement
from src.structure_baselines import all_canonical_7mers


ROOT = Path(__file__).resolve().parents[1]


def test_canonical_7mer_universe_has_8192_rc_classes():
    canonical = all_canonical_7mers()
    assert len(canonical) == 8192
    assert len(set(canonical)) == 8192
    assert all(canonical_rc(seq) == seq for seq in canonical)


def test_v0_4_scored_candidates_are_rc_class_units():
    scored = pd.read_parquet(ROOT / "data" / "processed" / "v0_4" / "v0_4_scored_candidates.parquet")
    counts = scored.groupby("protein_id")["canonical_7mer"].nunique()
    assert counts.eq(8192).all()
    assert not scored.duplicated(["protein_id", "canonical_7mer"]).any()
    assert scored["canonical_7mer"].map(canonical_rc).eq(scored["canonical_7mer"]).all()


def test_nampnn_predictions_do_not_duplicate_reverse_complement_classes():
    pred = pd.read_parquet(ROOT / "data" / "processed" / "v0_4" / "nampnn_predictions.parquet")
    assert set(pred["protein_id"]) == {"DBP35", "DBP48"}
    assert pred.groupby("protein_id")["canonical_7mer"].nunique().eq(8192).all()
    assert not pred.duplicated(["protein_id", "canonical_7mer"]).any()
    assert pred["canonical_7mer"].map(canonical_rc).eq(pred["canonical_7mer"]).all()
    assert pred["oriented_7mer"].map(reverse_complement).eq(pred["reverse_complement_7mer"]).all()


def test_split_helper_protects_reverse_complement_units():
    assert has_rc_split_leakage({"AAAAAAA"}, {"TTTTTTT"})
    assert not has_rc_split_leakage({"AAAAAAA"}, {"CCCCCCC"})


def test_v0_4_manifest_paths_are_project_relative():
    for rel_path in [
        ROOT / "metadata" / "v0_4" / "designed_dbp_structure_manifest.csv",
        ROOT / "metadata" / "v0_4" / "external_baseline_provenance.csv",
    ]:
        df = pd.read_csv(rel_path)
        for col in [c for c in df.columns if c.endswith("path") or c == "local_file" or c == "local_path"]:
            for value in df[col].dropna().astype(str):
                assert not Path(value).is_absolute()
                assert not PureWindowsPath(value).is_absolute()
                assert not value.startswith("E:")
