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

from src.utils import ensure_dir, project_root


ROOT = project_root()
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_FIGURES = ensure_dir(ROOT / "results" / "figures")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")


def save_gc_distribution(scored: pd.DataFrame) -> None:
    order = ["target", "single_mut", "double_mut", "gc_matched_random", "random_dna"]
    plot_df = scored[scored["candidate_type"].isin(order)].copy()
    plot_df["candidate_type"] = pd.Categorical(plot_df["candidate_type"], categories=order, ordered=True)
    summary = (
        plot_df.groupby("candidate_type")["gc_content"]
        .agg(["count", "mean", "std", "min", "max"])
        .reset_index()
    )
    summary.to_csv(RESULTS_TABLES / "gc_distribution_summary.csv", index=False)

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(data=plot_df, x="candidate_type", y="gc_content", inner="box", cut=0, linewidth=1, ax=ax, color="#d9e6f2")
    ax.set_xlabel("candidate type")
    ax.set_ylabel("GC content")
    ax.set_title("GC distribution across benchmark sets")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig2_gc_distribution.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    scored = pd.read_csv(PROCESSED_DIR / "benchmark_v0_1_scored.csv")
    save_gc_distribution(scored)
    order = ["target", "single_mut", "double_mut", "gc_matched_random", "random_dna"]
    plot_df = scored[scored["candidate_type"].isin(order)].copy()
    plot_df["candidate_type"] = pd.Categorical(plot_df["candidate_type"], categories=order, ordered=True)
    summary = plot_df.groupby("candidate_type")["sequence_proxy_score"].agg(["count", "mean", "std"]).reset_index()
    summary.to_csv(RESULTS_TABLES / "fig5_sequence_proxy_summary.csv", index=False)

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(data=plot_df, x="candidate_type", y="sequence_proxy_score", inner=None, cut=0, linewidth=1, ax=ax, color="#cfd8dc")
    sns.boxplot(data=plot_df, x="candidate_type", y="sequence_proxy_score", width=0.18, ax=ax, showfliers=False, color="#4c72b0")
    ax.set_xlabel("candidate type")
    ax.set_ylabel("sequence-only proxy score")
    ax.set_title("Target vs control sequence-only proxy distribution")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig5_sequence_proxy_distribution.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {RESULTS_FIGURES / 'fig5_specificity_proxy_distribution.png'}")


if __name__ == "__main__":
    main()
