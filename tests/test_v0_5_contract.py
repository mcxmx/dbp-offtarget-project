from pathlib import Path

import pandas as pd

from src.v0_5_contract import rank_order_unchanged_by_constant, stable_descending_order


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "metadata" / "v0_5"


def test_candidate_independent_constant_does_not_change_rank():
    scores = [0.4, 0.1, 0.4, -0.2]
    assert rank_order_unchanged_by_constant(scores, 17.5)
    assert stable_descending_order(scores) == stable_descending_order(
        [score - 17.5 for score in scores]
    )


def test_target_manifest_has_separate_target_concepts():
    manifest = pd.read_csv(METADATA / "designed_target_manifest_v0_5.csv")
    assert len(manifest) == 7
    required = {
        "intended_design_target",
        "experimental_assay_reference",
        "pbm_motif",
        "primary_target",
        "primary_target_source",
    }
    assert required.issubset(manifest.columns)
    dbp48 = manifest.loc[manifest["dbp_id"].eq("DBP48")].iloc[0]
    assert dbp48["intended_design_target"] == "CGCCCAAAGCCGCG"
    assert dbp48["experimental_assay_reference"] == "CGACACCTGACGCG"
    assert dbp48["pbm_motif"] == "CTGACG"
    assert dbp48["primary_target"] == dbp48["intended_design_target"]


def test_combined_components_reflect_known_target_and_protein_links():
    components = pd.read_csv(METADATA / "v0_5_split_components.csv")
    assert len(components) == 3
    members = set(components["members"])
    assert "DBP1|DBP3" in members
    assert "DBP35|DBP5|DBP6|DBP9" in members
    assert "DBP48" in members


def test_split_audit_reports_groups_and_legal_split_schemes():
    audit = pd.read_csv(METADATA / "v0_5_split_audit.csv")
    assert {"protein_cluster", "original_target_group", "assay_target_group", "motif_group"}.issubset(
        set(audit["record_type"])
    )
    schemes = audit.loc[audit["record_type"].eq("split_scheme")].set_index("group_id")
    assert schemes.loc["protein_cluster_loco", "legal_split"] == "yes"
    assert schemes.loc["protein_cluster_loco", "target_leakage_controlled"] == "no"
    assert schemes.loc["combined_component_loco", "target_leakage_controlled"] == "yes"
    assert schemes.loc["random_protein_7mer_row_split", "legal_split"] == "no"


def test_primary_split_never_separates_proteins_in_one_cluster():
    splits = pd.read_csv(METADATA / "v0_5_split_manifest.csv")
    primary = splits.loc[splits["split_name"].eq("protein_cluster_loco")]
    for (_, fold_id, cluster), group in primary.groupby(
        ["split_name", "fold_id", "protein_cluster"]
    ):
        assert group["partition"].nunique() == 1
    assert primary["row_level_random_split_allowed"].eq(False).all()
    assert primary["dna_unit"].eq("canonical_rc_equivalence_class").all()


def test_strict_split_keeps_combined_components_together():
    splits = pd.read_csv(METADATA / "v0_5_split_manifest.csv")
    strict = splits.loc[splits["split_name"].eq("combined_component_loco")]
    for (_, component), group in strict.groupby(["fold_id", "combined_component"]):
        assert group["partition"].nunique() == 1


def test_target_provenance_is_present_and_motif_is_not_primary():
    audit = pd.read_csv(METADATA / "target_provenance_audit_v0_5.csv")
    assert audit["dbp_id"].nunique() == 7
    assert set(audit["field"]) == {
        "intended_design_target",
        "experimental_assay_reference",
        "pbm_motif",
    }
    primary = audit.loc[audit["field"].eq("intended_design_target")]
    motif = audit.loc[audit["field"].eq("pbm_motif")]
    assert primary["independent_of_pbm"].eq(True).all()
    assert motif["independent_of_pbm"].eq(False).all()
