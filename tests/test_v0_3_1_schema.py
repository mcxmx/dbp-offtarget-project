from pathlib import Path, PureWindowsPath

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_v0_3_1_oriented_schema_uses_consensus_not_raw_name():
    oriented = pd.read_parquet(ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_oriented_v0_3_1.parquet")
    assert "experimental_score_raw" not in oriented.columns
    assert "experimental_escore_consensus" in oriented.columns
    assert "canonical_7mer" in oriented.columns
    assert oriented["oriented_7mer"].str.len().eq(7).all()
    assert oriented["canonical_7mer"].str.len().eq(7).all()


def test_rc_class_benchmark_counts():
    rc_class = pd.read_parquet(ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet")
    assert len(rc_class) == 7 * 8192
    assert rc_class.groupby("protein_id")["canonical_7mer"].nunique().eq(8192).all()
    assert rc_class["n_oriented_rows"].eq(2).all()
    assert rc_class["max_oriented_escore_abs_diff"].eq(0.0).all()


def test_portable_manifest_has_no_absolute_local_paths():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_3_1" / "gse237017_file_manifest_portable.csv")
    assert manifest["manifest_path_type"].eq("project_relative").all()
    for value in manifest["local_path"]:
        assert not Path(value).is_absolute()
        assert not PureWindowsPath(value).is_absolute()
        assert not str(value).startswith("E:")


def test_target_groups_and_clusters_are_present():
    groups = pd.read_csv(ROOT / "metadata" / "v0_3_1" / "designed_dbp_target_groups.csv")
    assert len(groups) == 7
    required = {"original_target_group", "assay_target_group", "motif_group", "protein_sequence_cluster"}
    assert required.issubset(groups.columns)
    assert not groups[list(required)].isna().any().any()
    assert groups["protein_sequence_cluster"].nunique() == 4


def test_disagreement_examples_are_not_total_counts():
    counts = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "all_disagreement_candidate_counts.csv")
    examples = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "top_disagreement_examples.csv")
    assert counts["n_disagreement"].sum() > len(examples)
    assert len(examples) == 7 * 20
