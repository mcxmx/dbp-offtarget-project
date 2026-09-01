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

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.utils import edit_distance, ensure_dir, kmer_jaccard, project_root


ROOT = project_root()
PROCESSED_DIR = ROOT / "data" / "processed" / "v0_3_1"
METADATA_DIR = ROOT / "metadata" / "v0_3_1"
TABLES_DIR = ensure_dir(ROOT / "results" / "v0_3_1" / "tables")
FIGURES_DIR = ensure_dir(ROOT / "results" / "v0_3_1" / "figures")
OLD_TABLES_DIR = ROOT / "results" / "v0_3" / "tables"


def spearman(x: pd.Series, y: pd.Series) -> float:
    return float(x.rank(method="average").corr(y.rank(method="average"), method="pearson"))


def windows(sequence: str, length: int) -> list[str]:
    if len(sequence) < length:
        return []
    return [sequence[i : i + length] for i in range(len(sequence) - length + 1)]


def best_window_similarity(candidate: str, motif: str) -> dict[str, float | str | bool]:
    motif_len = len(motif)
    motif_orientations = sorted({motif, reverse_complement(motif)})
    candidate_orientations = sorted({candidate, reverse_complement(candidate)})
    best_hamming = -1.0
    best_edit = -1.0
    best_window = ""
    best_reference = ""
    for candidate_orientation in candidate_orientations:
        for window in windows(candidate_orientation, motif_len):
            for reference in motif_orientations:
                mismatches = sum(a != b for a, b in zip(window, reference))
                hamming_similarity = 1.0 - mismatches / motif_len
                edit_similarity = 1.0 - edit_distance(window, reference) / motif_len
                if hamming_similarity > best_hamming:
                    best_hamming = hamming_similarity
                    best_edit = edit_similarity
                    best_window = window
                    best_reference = reference
    return {
        "hamming_similarity_to_paper_motif_rc_aware": best_hamming,
        "edit_similarity_to_paper_motif_rc_aware": best_edit,
        "best_candidate_window": best_window,
        "best_reference_orientation": best_reference,
        "motif_contained_rc_aware": any(
            reference in candidate_orientation
            for reference in motif_orientations
            for candidate_orientation in candidate_orientations
        ),
    }


def best_kmer_similarity(candidate: str, motif: str, k: int) -> float:
    values = []
    for candidate_orientation in sorted({candidate, reverse_complement(candidate)}):
        for motif_orientation in sorted({motif, reverse_complement(motif)}):
            values.append(kmer_jaccard(candidate_orientation, motif_orientation, k, rc_aware=True))
    return float(max(values)) if values else 0.0


def score_rc_classes() -> pd.DataFrame:
    benchmark = pd.read_parquet(PROCESSED_DIR / "designed_dbp_upbm_rc_class_v0_3_1.parquet")
    target_definitions = pd.read_csv(METADATA_DIR / "designed_dbp_target_definitions.csv")
    motif_by_protein = target_definitions.set_index("protein_id")["designed_binding_site_motif"].to_dict()
    rows = []
    for _, row in benchmark.iterrows():
        protein_id = row["protein_id"]
        candidate = row["canonical_7mer"]
        motif = motif_by_protein[protein_id]
        similarity = best_window_similarity(candidate, motif)
        rows.append(
            {
                "protein_id": protein_id,
                "canonical_7mer": candidate,
                "reverse_complement_7mer": reverse_complement(candidate),
                "motif_sequence": motif,
                "motif_canonical": canonical_rc(motif),
                "experimental_escore_consensus": row["experimental_escore_consensus"],
                **similarity,
                "kmer3_jaccard_to_paper_motif_rc_aware": best_kmer_similarity(candidate, motif, 3),
                "kmer4_jaccard_to_paper_motif_rc_aware": best_kmer_similarity(candidate, motif, 4),
                "sequence_unit": "protein_rc_equivalence_class",
                "notes": "RC-aware sequence-only proxy; not protein-conditioned.",
            }
        )
    return pd.DataFrame(rows)


def summarize_baseline(scored: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "hamming_similarity_to_paper_motif_rc_aware",
        "edit_similarity_to_paper_motif_rc_aware",
        "kmer3_jaccard_to_paper_motif_rc_aware",
        "kmer4_jaccard_to_paper_motif_rc_aware",
    ]
    rows = []
    for protein_id, group in scored.groupby("protein_id"):
        for metric in metrics:
            rows.append(
                {
                    "protein_id": protein_id,
                    "metric": metric,
                    "spearman": spearman(group[metric], group["experimental_escore_consensus"]),
                    "n_sequences": int(len(group)),
                    "sequence_unit": "RC equivalence class",
                    "reference_definition": "paper uPBM motif with candidate and motif reverse-complement equivalence",
                    "notes": "Per-protein sequence-only baseline; not protein-conditioned.",
                }
            )
    return pd.DataFrame(rows)


def compare_old_new(new_baseline: pd.DataFrame) -> pd.DataFrame:
    old = pd.read_csv(OLD_TABLES_DIR / "designed_dbp_sequence_baseline.csv")
    metric_map = {
        "hamming_similarity_to_target_7mer": "hamming_similarity_to_paper_motif_rc_aware",
        "edit_similarity_to_target_7mer": "edit_similarity_to_paper_motif_rc_aware",
        "kmer3_jaccard_to_target_7mer": "kmer3_jaccard_to_paper_motif_rc_aware",
        "kmer4_jaccard_to_target_7mer": "kmer4_jaccard_to_paper_motif_rc_aware",
    }
    rows = []
    for old_metric, new_metric in metric_map.items():
        old_sub = old[old["metric"] == old_metric].set_index("protein_id")
        new_sub = new_baseline[new_baseline["metric"] == new_metric].set_index("protein_id")
        for protein_id in sorted(set(old_sub.index) & set(new_sub.index)):
            old_value = float(old_sub.loc[protein_id, "spearman"])
            new_value = float(new_sub.loc[protein_id, "spearman"])
            rows.append(
                {
                    "protein_id": protein_id,
                    "metric": new_metric.replace("_similarity_to_paper_motif_rc_aware", "").replace("_jaccard_to_paper_motif_rc_aware", ""),
                    "old_metric": old_metric,
                    "new_metric": new_metric,
                    "old_spearman": old_value,
                    "rc_aware_spearman": new_value,
                    "delta": new_value - old_value,
                    "old_reference": "max similarity to intended-target-derived 7-mers; oriented v0.3 rows",
                    "new_reference": "paper motif; RC-class v0.3.1 rows",
                }
            )
    return pd.DataFrame(rows)


def plot_comparison(comparison: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    sns.boxplot(data=comparison, x="metric", y="delta", color="#d6e6f2", ax=ax)
    sns.stripplot(data=comparison, x="metric", y="delta", hue="protein_id", dodge=False, size=5, ax=ax)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("sequence-only metric")
    ax.set_ylabel("RC-aware Spearman minus v0.3 Spearman")
    ax.set_title("Effect of RC-aware motif baseline correction")
    ax.legend(title="protein", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_rc_aware_sequence_baseline_delta.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    scored = score_rc_classes()
    scored.to_parquet(PROCESSED_DIR / "designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1.parquet", index=False)
    scored.head(200).to_csv(PROCESSED_DIR / "designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1_preview.csv", index=False)
    baseline = summarize_baseline(scored)
    baseline.to_csv(TABLES_DIR / "designed_dbp_sequence_baseline_rc_aware.csv", index=False)
    comparison = compare_old_new(baseline)
    comparison.to_csv(TABLES_DIR / "designed_dbp_sequence_baseline_rc_aware_comparison.csv", index=False)
    plot_comparison(comparison)
    print(baseline.to_string(index=False))
    print(comparison.groupby("metric")[["old_spearman", "rc_aware_spearman", "delta"]].median().to_string())


if __name__ == "__main__":
    main()
