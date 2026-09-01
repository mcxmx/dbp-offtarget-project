from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.sequence_equivalence import canonical_rc
from src.utils import edit_distance, ensure_dir, project_root


ROOT = project_root()
METADATA_V03 = ROOT / "metadata" / "v0_3"
METADATA_DIR = ensure_dir(ROOT / "metadata" / "v0_3_1")
TABLES_DIR = ensure_dir(ROOT / "results" / "v0_3_1" / "tables")
FIGURES_DIR = ensure_dir(ROOT / "results" / "v0_3_1" / "figures")
CLUSTER_IDENTITY_THRESHOLD = 0.60


def simple_global_identity(seq1: str, seq2: str) -> float:
    # Transparent Levenshtein-style identity for short designed proteins.
    return 1.0 - edit_distance(seq1, seq2) / max(len(seq1), len(seq2), 1)


def connected_components(edges: dict[str, set[str]]) -> dict[str, str]:
    assignments = {}
    cluster_index = 1
    for node in sorted(edges):
        if node in assignments:
            continue
        stack = [node]
        members = []
        while stack:
            current = stack.pop()
            if current in assignments:
                continue
            assignments[current] = f"protein_cluster_{cluster_index}"
            members.append(current)
            stack.extend(sorted(edges[current] - set(assignments)))
        cluster_index += 1
    return assignments


def build_identity_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sequences = pd.read_csv(METADATA_V03 / "designed_dbp_sequences.csv")
    seq_by_protein = sequences.set_index("protein_id")["protein_sequence"].to_dict()
    proteins = sorted(seq_by_protein)
    matrix_rows = []
    long_rows = []
    edges = {protein: {protein} for protein in proteins}
    for protein_a in proteins:
        row = {"protein_id": protein_a}
        for protein_b in proteins:
            identity = simple_global_identity(seq_by_protein[protein_a], seq_by_protein[protein_b])
            row[protein_b] = identity
            long_rows.append(
                {
                    "protein_id_a": protein_a,
                    "protein_id_b": protein_b,
                    "sequence_identity": identity,
                    "identity_method": "1 - Levenshtein edit distance / max sequence length",
                    "cluster_threshold": CLUSTER_IDENTITY_THRESHOLD,
                }
            )
            if identity >= CLUSTER_IDENTITY_THRESHOLD:
                edges[protein_a].add(protein_b)
                edges[protein_b].add(protein_a)
        matrix_rows.append(row)
    assignments = connected_components(edges)
    cluster_rows = []
    for protein in proteins:
        cluster_rows.append(
            {
                "protein_id": protein,
                "protein_sequence_cluster": assignments[protein],
                "cluster_identity_threshold": CLUSTER_IDENTITY_THRESHOLD,
                "cluster_method": "connected components using simple_global_identity >= threshold",
            }
        )
    return pd.DataFrame(matrix_rows), pd.DataFrame(long_rows), pd.DataFrame(cluster_rows)


def build_target_groups(clusters: pd.DataFrame) -> pd.DataFrame:
    target_definitions = pd.read_csv(METADATA_DIR / "designed_dbp_target_definitions.csv")
    cluster_by_protein = clusters.set_index("protein_id")["protein_sequence_cluster"].to_dict()
    rows = []
    for _, row in target_definitions.sort_values("protein_id").iterrows():
        motif = row["designed_binding_site_motif"]
        rows.append(
            {
                "protein_id": row["protein_id"],
                "original_target_group": f"original_target_{row['original_design_target_id']}",
                "assay_target_group": f"assay_target_{row['experimental_assay_target_id']}",
                "motif_group": f"motif_{canonical_rc(motif)}",
                "designed_binding_site_motif": motif,
                "protein_sequence_cluster": cluster_by_protein[row["protein_id"]],
                "notes": "Future splits should not treat shared target/motif/sequence-cluster proteins as independent strongest generalization tests.",
            }
        )
    return pd.DataFrame(rows)


def plot_identity(matrix: pd.DataFrame) -> None:
    heatmap = matrix.set_index("protein_id")
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    sns.heatmap(heatmap, annot=True, fmt=".2f", cmap="viridis", vmin=0, vmax=1, cbar_kws={"label": "sequence identity"}, ax=ax)
    ax.set_xlabel("designed DBP")
    ax.set_ylabel("designed DBP")
    ax.set_title("Designed DBP protein sequence similarity")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_designed_dbp_sequence_similarity.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    matrix, long, clusters = build_identity_tables()
    matrix.to_csv(TABLES_DIR / "designed_dbp_sequence_identity_matrix.csv", index=False)
    long.to_csv(TABLES_DIR / "designed_dbp_pairwise_sequence_identity.csv", index=False)
    clusters.to_csv(TABLES_DIR / "designed_dbp_sequence_clusters.csv", index=False)
    target_groups = build_target_groups(clusters)
    target_groups.to_csv(METADATA_DIR / "designed_dbp_target_groups.csv", index=False)
    plot_identity(matrix)
    print(matrix.to_string(index=False))
    print(clusters.to_string(index=False))
    print(target_groups.to_string(index=False))


if __name__ == "__main__":
    main()
