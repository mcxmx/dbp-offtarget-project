from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import compute_sequence_metrics, ensure_dir, project_root


ROOT = project_root()
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_FIGURES = ensure_dir(ROOT / "results" / "figures")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")


METRICS = [
    "sequence_identity",
    "edit_similarity",
    "hamming_similarity",
    "kmer3_jaccard",
    "kmer4_jaccard",
    "rc_kmer4_jaccard",
]


def add_consensus_similarity(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source_id, group in df.groupby("source_id"):
        consensus_rows = group[group["sequence_class"] == "consensus"]
        if consensus_rows.empty:
            continue
        consensus = consensus_rows.iloc[0]["dna_sequence"]
        for _, row in group.iterrows():
            metrics = compute_sequence_metrics(consensus, row["dna_sequence"])
            length = max(len(consensus), 1)
            rows.append(
                {
                    **row.to_dict(),
                    "consensus_sequence": consensus,
                    "sequence_identity": metrics["sequence_identity"],
                    "edit_similarity": 1.0 - metrics["edit_distance"] / length,
                    "hamming_similarity": 1.0 - metrics["hamming_distance"] / length,
                    "kmer3_jaccard": metrics["kmer3_jaccard"],
                    "kmer4_jaccard": metrics["kmer4_jaccard"],
                    "rc_kmer4_jaccard": metrics["rc_kmer4_jaccard"],
                }
            )
    return pd.DataFrame(rows)


def spearman_summary(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (source_id, protein_name), group in scored.groupby(["source_id", "protein_name"]):
        for metric in METRICS:
            metric_rank = group[metric].rank(method="average")
            score_rank = group["experimental_score"].rank(method="average")
            rows.append(
                {
                    "source_id": source_id,
                    "protein_name": protein_name,
                    "metric": metric,
                    "n_sequences": len(group),
                    "spearman_with_pfm_derived_score": metric_rank.corr(score_rank, method="pearson"),
                    "score_type": group["score_type"].iloc[0],
                    "notes": "Correlation uses sequence similarity to the JASPAR consensus versus PFM-derived PWM log2-odds score.",
                }
            )
    return pd.DataFrame(rows)


def plot(summary: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    plot_df = summary.copy()
    plot_df["metric"] = pd.Categorical(plot_df["metric"], categories=METRICS, ordered=True)
    sns.barplot(
        data=plot_df,
        x="metric",
        y="spearman_with_pfm_derived_score",
        hue="protein_name",
        ax=ax,
    )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("sequence-only metric to consensus")
    ax.set_ylabel("Spearman correlation")
    ax.set_title("PFM-derived score vs sequence-only proxies")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="protein", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig_experimental_specificity_proxy_baseline.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    df = pd.read_csv(PROCESSED_DIR / "experimental_specificity_small.csv")
    scored = add_consensus_similarity(df)
    scored.to_csv(PROCESSED_DIR / "experimental_specificity_small_with_sequence_proxy.csv", index=False)
    summary = spearman_summary(scored)
    summary.to_csv(RESULTS_TABLES / "experimental_specificity_proxy_spearman.csv", index=False)
    plot(summary)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
