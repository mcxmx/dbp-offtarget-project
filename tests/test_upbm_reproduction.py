from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_paper_percentile_reproduction_within_tolerance():
    reproduction = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "paper_percentile_reproduction.csv")
    assert len(reproduction) == 7
    assert reproduction["reproduction_status"].eq("PASS").all()
    assert (reproduction["absolute_difference"] <= reproduction["qc_tolerance_percentile_points"]).all()
    assert set(reproduction["protein_id"]) == {"DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"}


def test_dbp48_target_definitions_are_separated():
    targets = pd.read_csv(ROOT / "metadata" / "v0_3_1" / "designed_dbp_target_definitions.csv")
    dbp48 = targets.set_index("protein_id").loc["DBP48"]
    assert dbp48["original_design_target_id"] == "I_b"
    assert dbp48["experimental_assay_target_id"] == "C"
    assert dbp48["original_design_target"] != dbp48["experimental_assay_target"]
    assert dbp48["designed_binding_site_motif"] == "CTGACG"


def test_paper_source_data_manifest_has_hashes():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_3_1" / "paper_source_data_manifest.csv")
    assert len(manifest) == 2
    assert manifest["sha256"].str.fullmatch("[0-9a-f]{64}").all()
    assert manifest["download_status"].eq("downloaded").all()
