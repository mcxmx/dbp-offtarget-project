from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed" / "v0_3")
TABLES_DIR = ensure_dir(ROOT / "results" / "v0_3" / "tables")
FIGURES_DIR = ensure_dir(ROOT / "results" / "v0_3" / "figures")
BASES = "ACGT"


def fig_escore_distribution(benchmark: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.violinplot(
        data=benchmark,
        x="protein_id",
        y="experimental_score_primary",
        inner=None,
        cut=0,
        linewidth=1,
        color="#d9e6f2",
        ax=ax,
    )
    sns.boxplot(
        data=benchmark,
        x="protein_id",
        y="experimental_score_primary",
        width=0.18,
        showfliers=False,
        color="#4c72b0",
        ax=ax,
    )
    ax.set_xlabel("designed DBP")
    ax.set_ylabel("experimental PBM E-score")
    ax.set_title("uPBM E-score distributions by designed DBP")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_v0_3_3_escore_distributions_by_dbp.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def top_kmer_base_frequency(benchmark: pd.DataFrame) -> pd.DataFrame:
    rows = []
    top_rows = []
    for protein_id, group in benchmark.groupby("protein_id"):
        top = group[group["experimental_percentile"] >= 0.99].copy()
        if top.empty:
            top = group.nlargest(max(1, int(np.ceil(0.01 * len(group)))), "experimental_score_primary").copy()
        top_rows.append(top[["protein_id", "dna_7mer", "experimental_score_primary", "experimental_percentile"]])
        for position in range(7):
            bases = top["dna_7mer"].str[position]
            counts = bases.value_counts(normalize=True)
            for base in BASES:
                rows.append(
                    {
                        "protein_id": protein_id,
                        "position": position + 1,
                        "base": base,
                        "frequency_in_top_1_percent": float(counts.get(base, 0.0)),
                        "n_top_7mers": int(len(top)),
                    }
                )
    pd.concat(top_rows, ignore_index=True).to_csv(TABLES_DIR / "top_1_percent_7mers_by_protein.csv", index=False)
    return pd.DataFrame(rows)


def fig_top_kmer_landscape(base_freq: pd.DataFrame) -> None:
    proteins = sorted(base_freq["protein_id"].unique())
    sns.set_theme(style="white", context="paper")
    fig, axes = plt.subplots(len(proteins), 1, figsize=(9, 1.45 * len(proteins)), sharex=True)
    if len(proteins) == 1:
        axes = [axes]
    for ax, protein_id in zip(axes, proteins):
        sub = base_freq[base_freq["protein_id"] == protein_id]
        matrix = sub.pivot(index="base", columns="position", values="frequency_in_top_1_percent").loc[list(BASES)]
        sns.heatmap(matrix, ax=ax, cmap="viridis", vmin=0, vmax=1, cbar=protein_id == proteins[-1], cbar_kws={"label": "base frequency"})
        ax.set_ylabel(protein_id)
        ax.set_xlabel("")
    axes[-1].set_xlabel("7-mer position")
    fig.suptitle("Top 1% experimental 7-mer base preferences", y=1.0)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_v0_3_4_top_7mer_base_preferences.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_target_rank(target_summary: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    sns.barplot(data=target_summary, x="protein_id", y="best_target_percentile", ax=axes[0], color="#55a868")
    axes[0].set_xlabel("designed DBP")
    axes[0].set_ylabel("best target-derived 7-mer percentile")
    axes[0].set_title("Best intended-target 7-mer rank")
    axes[0].set_ylim(0, 1)
    sns.barplot(data=target_summary, x="protein_id", y="mean_target_7mer_escore", ax=axes[1], color="#c44e52")
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xlabel("designed DBP")
    axes[1].set_ylabel("mean target-derived 7-mer E-score")
    axes[1].set_title("Target-derived 7-mer mean E-score")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_v0_3_5_intended_target_7mer_rank_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_sequence_similarity_baseline(baseline: pd.DataFrame) -> None:
    matrix = baseline.pivot(index="protein_id", columns="metric", values="spearman")
    metric_order = [
        "hamming_similarity_to_target_7mer",
        "edit_similarity_to_target_7mer",
        "kmer3_jaccard_to_target_7mer",
        "kmer4_jaccard_to_target_7mer",
    ]
    matrix = matrix[metric_order]
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap="vlag", center=0, vmin=-0.3, vmax=0.3, cbar_kws={"label": "Spearman"})
    ax.set_xlabel("sequence-only metric to target-derived 7-mers")
    ax.set_ylabel("designed DBP")
    ax.set_title("Sequence similarity vs experimental PBM E-score")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_v0_3_6_sequence_similarity_vs_escore.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def failure_candidates(benchmark: pd.DataFrame, baseline_scored: pd.DataFrame) -> pd.DataFrame:
    scored = baseline_scored.merge(
        benchmark[["protein_id", "dna_7mer", "experimental_percentile", "experimental_rank"]],
        on=["protein_id", "dna_7mer"],
        how="left",
    )
    rows = []
    for protein_id, group in scored.groupby("protein_id"):
        median_sim = group["hamming_similarity_to_target_7mer"].median()
        candidates = group[
            (group["experimental_percentile"] >= 0.99)
            & (group["hamming_similarity_to_target_7mer"] <= median_sim)
        ].copy()
        candidates = candidates.sort_values("experimental_score_primary", ascending=False).head(20)
        candidates["protein_median_hamming_similarity_to_target_7mer"] = median_sim
        candidates["selection_rule"] = "top 1% PBM E-score and hamming similarity to target-derived 7-mers at or below protein median"
        rows.append(candidates)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> None:
    benchmark = pd.read_parquet(PROCESSED_DIR / "designed_dbp_upbm_v0_3.parquet")
    baseline = pd.read_csv(TABLES_DIR / "designed_dbp_sequence_baseline.csv")
    target_summary = pd.read_csv(TABLES_DIR / "target_rank_summary.csv")
    baseline_scored = pd.read_parquet(PROCESSED_DIR / "designed_dbp_sequence_baseline_scored_v0_3.parquet")

    fig_escore_distribution(benchmark)
    base_freq = top_kmer_base_frequency(benchmark)
    base_freq.to_csv(TABLES_DIR / "top_1_percent_7mer_base_frequencies.csv", index=False)
    fig_top_kmer_landscape(base_freq)
    fig_target_rank(target_summary)
    fig_sequence_similarity_baseline(baseline)
    failures = failure_candidates(benchmark, baseline_scored)
    failures.to_csv(TABLES_DIR / "sequence_similarity_failure_candidates.csv", index=False)
    print(f"figures written to {FIGURES_DIR}")
    print(f"sequence-similarity failure candidates: {len(failures)}")
    print(baseline.groupby("metric")["spearman"].agg(["min", "median", "max"]).to_string())


if __name__ == "__main__":
    main()
