from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "metadata" / "v1_0_structure"
OUT = ROOT / "results" / "v1_0_structure"


def test_structure_inventory_has_all_seven_exact_design_sequences():
    inv = pd.read_csv(META / "structure_inventory.csv")
    assert list(inv["protein_id"]) == ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]
    assert inv["structure_parse_status"].eq("parsed").all()
    assert inv["exact_pbm_construct_in_structure"].eq(True).all()
    assert inv["protein_chain_id"].eq("A").all()
    assert inv["dna_chain_ids"].eq("B").all()


def test_coverage_audit_does_not_call_old_coverage_biological():
    audit = pd.read_csv(OUT / "v0_8_coverage_failure_audit.csv")
    assert len(audit) == 14
    assert audit["biological_limitation_inferred"].eq(False).all()
    assert set(audit["coverage_failure_class"]) == {"previous_pipeline_resource_audit_gap"}


def test_nampnn_full_landscape_is_protein_level_complete():
    metrics = pd.read_csv(OUT / "nampnn_metrics.csv")
    assert metrics["status"].eq("evaluated").all()
    assert metrics["n_candidates"].eq(8192).all()
    pred = pd.read_csv(OUT / "nampnn_results.csv")
    assert pred["protein_id"].nunique() == 7
    assert pred.groupby("protein_id")["canonical_7mer"].nunique().eq(8192).all()


def test_rossetta_and_af3_subsets_are_frozen_before_execution():
    ros = pd.read_csv(OUT / "rossetta_subset_results.csv")
    af3 = pd.read_csv(OUT / "af3_contact_subset_results.csv")
    assert set(ros["protein_id"]) == {"DBP35", "DBP48"}
    assert set(af3["protein_id"]) == {"DBP35", "DBP48"}
    assert ros["selection_rule"].str.contains("before").all()
    assert af3["selection_rule"].str.contains("before").all()
    assert ros["status"].str.startswith("not_run").all()
    assert af3["status"].str.startswith("not_run").all()


def test_v0_9_external_firewall_remains_locked():
    manifest = (ROOT / "metadata" / "v0_9_external" / "FROZEN_VALIDATION_MANIFEST.yaml").read_text()
    assert "label_access: false" in manifest
