from pathlib import Path, PureWindowsPath

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _assert_no_absolute_paths(df: pd.DataFrame) -> None:
    for col in df.columns:
        if df[col].dtype != object:
            continue
        for value in df[col].dropna().astype(str):
            assert not Path(value).is_absolute()
            assert not PureWindowsPath(value).is_absolute()
            assert not value.startswith("E:")


def test_external_baseline_provenance_records_versions_and_status():
    provenance = pd.read_csv(ROOT / "metadata" / "v0_4" / "external_baseline_provenance.csv")
    assert {"DeepPBS", "NA-MPNN", "SimpleProteinConditionalBaseline"}.issubset(set(provenance["baseline"]))
    repo_rows = provenance[provenance["baseline"].isin(["DeepPBS", "NA-MPNN"])]
    assert repo_rows["commit_hash"].str.fullmatch("[0-9a-f]{40}").all()
    assert provenance["inference_status_v0_4"].notna().all()
    assert "not_evaluable_in_current_environment" in set(provenance["inference_status_v0_4"])
    assert "ran_for_dbp35_and_dbp48_only" in set(provenance["inference_status_v0_4"])
    _assert_no_absolute_paths(provenance)


def test_overlap_audit_flags_dbp48_nampnn_not_zero_shot():
    overlap = pd.read_csv(ROOT / "metadata" / "v0_4" / "baseline_data_overlap_audit.csv")
    assert not overlap.empty
    dbp48 = overlap[(overlap["protein_id"] == "DBP48") & (overlap["baseline"] == "NA-MPNN")].iloc[0]
    assert bool(dbp48["exact_protein_seen"])
    assert bool(dbp48["structure_seen"])
    assert dbp48["risk_level"] == "high_not_zero_shot"
    assert "8tac" in dbp48["evidence"].lower()


def test_structure_manifest_records_evaluable_and_missing_cases():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_4" / "designed_dbp_structure_manifest.csv")
    assert set(manifest["protein_id"]) == {"DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"}
    assert manifest.loc[manifest["protein_id"] == "DBP35", "local_file"].iloc[0] == "external/dbp_design/2b_design_mpnn/DBP035.pdb"
    assert manifest.loc[manifest["protein_id"] == "DBP48", "pdb_id"].iloc[0] == "8TAC"
    assert (manifest["structure_confidence"] == "none").sum() == 5
    _assert_no_absolute_paths(manifest)


def test_v0_4_prediction_paths_are_portable():
    pred = pd.read_parquet(ROOT / "data" / "processed" / "v0_4" / "nampnn_predictions.parquet")
    assert "source_npz" in pred.columns
    _assert_no_absolute_paths(pred[["source_npz"]])
