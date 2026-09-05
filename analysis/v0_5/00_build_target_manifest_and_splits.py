from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.v0_5_contract import build_components
from src.utils import ensure_dir, project_root


ROOT = project_root()
V03_METADATA = ROOT / "metadata" / "v0_3"
V031_METADATA = ROOT / "metadata" / "v0_3_1"
V031_TABLES = ROOT / "results" / "v0_3_1" / "tables"
OUT_METADATA = ensure_dir(ROOT / "metadata" / "v0_5")

DBP_ORDER = ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]


def load_inputs() -> pd.DataFrame:
    sequences = pd.read_csv(V03_METADATA / "designed_dbp_sequences.csv")
    targets = pd.read_csv(V031_METADATA / "designed_dbp_target_definitions.csv")
    groups = pd.read_csv(V031_METADATA / "designed_dbp_target_groups.csv")
    clusters = pd.read_csv(V031_TABLES / "designed_dbp_sequence_clusters.csv")

    required = {
        "protein_id",
        "original_design_target",
        "original_design_target_source",
        "experimental_assay_target",
        "experimental_assay_target_source",
        "designed_binding_site_motif",
        "motif_source",
        "source_url",
        "paper_doi",
        "retrieval_date",
    }
    missing = required - set(targets.columns)
    if missing:
        raise ValueError(f"Target definition columns missing: {sorted(missing)}")

    merged = (
        sequences[["protein_id", "protein_sequence"]]
        .merge(
            targets[
                [
                    "protein_id",
                    "original_design_target_id",
                    "original_design_target",
                    "original_design_target_source",
                    "experimental_assay_target_id",
                    "experimental_assay_target",
                    "experimental_assay_target_source",
                    "designed_binding_site_motif",
                    "designed_binding_site_motif_length",
                    "motif_source",
                    "source_url",
                    "paper_doi",
                    "retrieval_date",
                    "confidence",
                ]
            ],
            on="protein_id",
            validate="one_to_one",
        )
        .merge(
            groups[
                [
                    "protein_id",
                    "original_target_group",
                    "assay_target_group",
                    "motif_group",
                    "protein_sequence_cluster",
                ]
            ],
            on="protein_id",
            validate="one_to_one",
        )
        .merge(
            clusters[["protein_id", "cluster_identity_threshold", "cluster_method"]],
            on="protein_id",
            validate="one_to_one",
        )
    )
    merged["protein_sequence_cluster"] = merged["protein_sequence_cluster"].astype(str)
    merged = merged.set_index("protein_id").loc[DBP_ORDER].reset_index()
    if merged["protein_id"].duplicated().any() or len(merged) != len(DBP_ORDER):
        raise ValueError("Expected exactly one row for each of the seven designed DBPs")
    return merged


def build_target_manifest(inputs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in inputs.to_dict("records"):
        original = row["original_design_target"]
        assay = row["experimental_assay_target"]
        motif = row["designed_binding_site_motif"]
        target_group = "|".join(
            [
                row["original_target_group"],
                row["assay_target_group"],
                row["motif_group"],
            ]
        )
        primary_source = row["original_design_target_source"]
        notes = (
            "primary_target is the independently reported original design target; "
            "it is not inferred from PBM scores or top PBM k-mers. "
            "The assay reference and PBM-derived motif remain separate fields."
        )
        if row["protein_id"] == "DBP48":
            notes += (
                " DBP48 preserves original design target I_b, assay/reference target C, "
                "and PBM motif CTGACG as distinct concepts."
            )
        rows.append(
            {
                "dbp_id": row["protein_id"],
                "protein_sequence": row["protein_sequence"],
                "protein_cluster": row["protein_sequence_cluster"],
                "intended_design_target": original,
                "intended_target_source": primary_source,
                "experimental_assay_reference": assay,
                "assay_reference_source": row["experimental_assay_target_source"],
                "pbm_motif": motif,
                "target_group": target_group,
                "original_target_group": row["original_target_group"],
                "assay_target_group": row["assay_target_group"],
                "motif_group": row["motif_group"],
                "primary_target": original,
                "primary_target_source": primary_source,
                "notes": notes,
                "source_url": row["source_url"],
                "paper_doi": row["paper_doi"],
                "retrieval_date": row["retrieval_date"],
                "target_confidence": row["confidence"],
            }
        )
    return pd.DataFrame(rows)


def build_provenance_audit(inputs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    fields = [
        (
            "intended_design_target",
            "original_design_target",
            "original_design_target_source",
            True,
            "independent design metadata",
        ),
        (
            "experimental_assay_reference",
            "experimental_assay_target",
            "experimental_assay_target_source",
            True,
            "independent assay/reference metadata",
        ),
        (
            "pbm_motif",
            "designed_binding_site_motif",
            "motif_source",
            False,
            "PBM evaluation motif; not used as primary target",
        ),
    ]
    for row in inputs.to_dict("records"):
        for field, source_column, source_column_for_provenance, independent, note in fields:
            rows.append(
                {
                    "dbp_id": row["protein_id"],
                    "field": field,
                    "value": row[source_column],
                    "source_file": "metadata/v0_3_1/designed_dbp_target_definitions.csv",
                    "source_column_or_rule": source_column_for_provenance,
                    "provenance_status": "resolved",
                    "independent_of_pbm": independent,
                    "source_url": row["source_url"],
                    "paper_doi": row["paper_doi"],
                    "retrieval_date": row["retrieval_date"],
                    "notes": note,
                }
            )
    return pd.DataFrame(rows)


def build_leakage_graph(
    manifest: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    edges = []
    for left, right in combinations(manifest["dbp_id"], 2):
        left_row = manifest.loc[manifest["dbp_id"].eq(left)].iloc[0]
        right_row = manifest.loc[manifest["dbp_id"].eq(right)].iloc[0]
        same_cluster = left_row["protein_cluster"] == right_row["protein_cluster"]
        same_original = left_row["original_target_group"] == right_row["original_target_group"]
        same_assay = left_row["assay_target_group"] == right_row["assay_target_group"]
        same_motif = left_row["motif_group"] == right_row["motif_group"]
        reasons = []
        if same_cluster:
            reasons.append("protein_cluster")
        if same_original:
            reasons.append("original_target_group")
        if same_assay:
            reasons.append("assay_target_group")
        if same_motif:
            reasons.append("motif_group")
        if reasons:
            edges.append(
                {
                    "dbp_id_a": left,
                    "dbp_id_b": right,
                    "same_protein_cluster": same_cluster,
                    "same_original_target_group": same_original,
                    "same_assay_target_group": same_assay,
                    "same_motif_group": same_motif,
                    "relationship_basis": "|".join(reasons),
                    "target_leakage_risk": "yes" if any([same_original, same_assay, same_motif]) else "no",
                }
            )

    edge_df = pd.DataFrame(edges)
    component_map = build_components(
        manifest["dbp_id"],
        edge_df[["dbp_id_a", "dbp_id_b"]].itertuples(index=False, name=None),
    )
    components = sorted(component_map.values(), key=lambda members: members[0])
    component_id_by_protein = {}
    component_rows = []
    for index, members in enumerate(components, start=1):
        component_id = f"combined_component_{index}"
        for protein in members:
            component_id_by_protein[protein] = component_id
        subset = manifest[manifest["dbp_id"].isin(members)]
        component_rows.append(
            {
                "component_id": component_id,
                "members": "|".join(members),
                "n_proteins": len(members),
                "protein_clusters": "|".join(sorted(subset["protein_cluster"].unique())),
                "original_target_groups": "|".join(sorted(subset["original_target_group"].unique())),
                "assay_target_groups": "|".join(sorted(subset["assay_target_group"].unique())),
                "motif_groups": "|".join(sorted(subset["motif_group"].unique())),
                "legal_for_primary_cluster_split": True,
                "legal_for_strict_component_split": True,
                "notes": "All proteins in a component must remain in the same future split.",
            }
        )
    edge_df["combined_component_a"] = edge_df["dbp_id_a"].map(component_id_by_protein)
    edge_df["combined_component_b"] = edge_df["dbp_id_b"].map(component_id_by_protein)
    edge_df["combined_component"] = edge_df["combined_component_a"]
    edge_df = edge_df.drop(columns=["combined_component_a", "combined_component_b"])
    component_df = pd.DataFrame(component_rows)

    audit_rows = []
    for column, record_type in [
        ("protein_cluster", "protein_cluster"),
        ("original_target_group", "original_target_group"),
        ("assay_target_group", "assay_target_group"),
        ("motif_group", "motif_group"),
    ]:
        for group_id, subset in manifest.groupby(column, sort=True):
            audit_rows.append(
                {
                    "record_type": record_type,
                    "group_id": group_id,
                    "members": "|".join(sorted(subset["dbp_id"])),
                    "n_members": len(subset),
                    "protein_clusters": "|".join(sorted(subset["protein_cluster"].unique())),
                    "original_target_groups": "|".join(sorted(subset["original_target_group"].unique())),
                    "assay_target_groups": "|".join(sorted(subset["assay_target_group"].unique())),
                    "motif_groups": "|".join(sorted(subset["motif_group"].unique())),
                    "legal_split": "primary_cluster_loco_and_strict_component_loco",
                    "target_leakage_controlled": "strict_component_loco_only",
                    "evidence": "existing v0.3.1 curated metadata",
                    "notes": "Group members must remain together for any split that claims independence at this level.",
                }
            )
    for row in component_df.to_dict("records"):
        audit_rows.append(
            {
                "record_type": "combined_component",
                "group_id": row["component_id"],
                "members": row["members"],
                "n_members": row["n_proteins"],
                "protein_clusters": row["protein_clusters"],
                "original_target_groups": row["original_target_groups"],
                "assay_target_groups": row["assay_target_groups"],
                "motif_groups": row["motif_groups"],
                "legal_split": "strict_component_loco",
                "target_leakage_controlled": "yes",
                "evidence": "connected components of protein/target/motif leakage graph",
                "notes": row["notes"],
            }
        )
    split_rows = [
        (
            "protein_cluster_loco",
            "yes",
            "no",
            "primary protein-level split; target/motif groups can still cross folds",
        ),
        (
            "combined_component_loco",
            "yes",
            "yes",
            "strict sensitivity split; all known protein/target/motif edges stay within folds",
        ),
        (
            "random_protein_7mer_row_split",
            "no",
            "no",
            "prohibited because it leaks protein identity and can split related target/motif groups",
        ),
    ]
    for scheme, legal, target_controlled, notes in split_rows:
        audit_rows.append(
            {
                "record_type": "split_scheme",
                "group_id": scheme,
                "members": "",
                "n_members": pd.NA,
                "protein_clusters": "",
                "original_target_groups": "",
                "assay_target_groups": "",
                "motif_groups": "",
                "legal_split": legal,
                "target_leakage_controlled": target_controlled,
                "evidence": "v0.5 split contract",
                "notes": notes,
            }
        )
    audit_df = pd.DataFrame(audit_rows)
    return edge_df, component_df, audit_df


def build_split_manifest(manifest: pd.DataFrame, components: pd.DataFrame) -> pd.DataFrame:
    rows = []
    schemes = [
        (
            "protein_cluster_loco",
            "primary",
            "protein_cluster",
            sorted(manifest["protein_cluster"].unique()),
        ),
        (
            "combined_component_loco",
            "strict_sensitivity",
            "combined_component",
            sorted(components["component_id"].unique()),
        ),
    ]
    component_by_protein = {}
    for component in components.to_dict("records"):
        for protein in component["members"].split("|"):
            component_by_protein[protein] = component["component_id"]

    for split_name, split_role, split_key, held_out_values in schemes:
        for fold_index, held_out in enumerate(held_out_values, start=1):
            fold_id = f"{split_name}_fold_{fold_index}"
            for row in manifest.to_dict("records"):
                protein = row["dbp_id"]
                protein_value = (
                    row["protein_cluster"]
                    if split_key == "protein_cluster"
                    else component_by_protein[protein]
                )
                partition = "test" if protein_value == held_out else "train"
                rows.append(
                    {
                        "split_name": split_name,
                        "split_role": split_role,
                        "fold_id": fold_id,
                        "held_out_group": held_out,
                        "dbp_id": protein,
                        "partition": partition,
                        "protein_cluster": row["protein_cluster"],
                        "combined_component": component_by_protein[protein],
                        "original_target_group": row["original_target_group"],
                        "assay_target_group": row["assay_target_group"],
                        "motif_group": row["motif_group"],
                        "dna_unit": "canonical_rc_equivalence_class",
                        "row_level_random_split_allowed": False,
                        "validation_partition": "not_defined_in_v0_5_contract",
                        "notes": (
                            "Future model training must split at protein level; DNA units remain "
                            "grouped by canonical reverse-complement class."
                        ),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    inputs = load_inputs()
    manifest = build_target_manifest(inputs)
    provenance = build_provenance_audit(inputs)
    edge_df, components, audit_df = build_leakage_graph(manifest)
    split_manifest = build_split_manifest(manifest, components)

    manifest.to_csv(OUT_METADATA / "designed_target_manifest_v0_5.csv", index=False)
    provenance.to_csv(OUT_METADATA / "target_provenance_audit_v0_5.csv", index=False)
    audit_df.to_csv(OUT_METADATA / "v0_5_split_audit.csv", index=False)
    edge_df.to_csv(OUT_METADATA / "v0_5_split_edges.csv", index=False)
    components.to_csv(OUT_METADATA / "v0_5_split_components.csv", index=False)
    split_manifest.to_csv(OUT_METADATA / "v0_5_split_manifest.csv", index=False)

    print(f"target manifest: {len(manifest)} proteins")
    print(f"split edges: {len(edge_df)}")
    print(f"combined components: {len(components)}")
    print(components[["component_id", "members"]].to_string(index=False))
    print(
        split_manifest.groupby(["split_name", "fold_id", "partition"])["dbp_id"]
        .nunique()
        .to_string()
    )


if __name__ == "__main__":
    main()
