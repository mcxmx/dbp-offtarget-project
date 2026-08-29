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
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import compute_sequence_metrics, ensure_dir, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_FIGURES = ensure_dir(ROOT / "results" / "figures" / "v0_2")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")


def resolve_pair(pairs: pd.DataFrame, key: str) -> pd.Series:
    exact = pairs.loc[pairs["pair_id"] == key]
    if not exact.empty:
        return exact.iloc[0]
    by_pdb = pairs.loc[pairs["pdb_id"] == key]
    if not by_pdb.empty:
        return by_pdb.iloc[0]
    prefix = pairs.loc[pairs["pair_id"].astype(str).str.startswith(f"{key}_")]
    if not prefix.empty:
        return prefix.iloc[0]
    raise ValueError(f"Pair not found: {key}")


def fig1_dataset_overview(benchmark: pd.DataFrame, pairs: pd.DataFrame) -> None:
    counts = (
        benchmark.groupby(["pair_id", "candidate_type"])
        .size()
        .reset_index(name="count")
        .pivot(index="pair_id", columns="candidate_type", values="count")
        .fillna(0)
    )
    order = pairs.sort_values("dna_length", ascending=False)["pair_id"].tolist()
    candidate_order = ["target", "single_mut", "double_mut", "gc_matched_random", "random_dna"]
    for col in candidate_order:
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[candidate_order].reindex(order).fillna(0)
    counts.to_csv(RESULTS_TABLES / "dataset_overview_counts_v0_2.csv")

    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [3.1, 1.2]})
    sns.heatmap(
        np.log10(counts + 1),
        ax=axes[0],
        cmap="viridis",
        cbar_kws={"label": "log10(count + 1)"},
        linewidths=0.3,
        linecolor="white",
    )
    axes[0].set_xlabel("candidate type")
    axes[0].set_ylabel("curated pair_id")
    axes[0].set_title("v0.2 curated benchmark overview")

    lengths = pairs.set_index("pair_id").loc[order, "dna_length"]
    axes[1].barh(range(len(order)), lengths.values, color="#4c72b0")
    axes[1].set_yticks(range(len(order)))
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("target DNA length")
    axes[1].set_title("Target length")
    axes[1].invert_yaxis()
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig1_dataset_overview_v0_2.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig2_gc_distribution(scored: pd.DataFrame) -> None:
    order = ["target", "single_mut", "double_mut", "gc_matched_random", "random_dna"]
    plot_df = scored[scored["candidate_type"].isin(order)].copy()
    plot_df["candidate_type"] = pd.Categorical(plot_df["candidate_type"], categories=order, ordered=True)
    summary = (
        plot_df.groupby("candidate_type", observed=False)["gc_content"]
        .agg(["count", "mean", "std", "min", "max"])
        .reset_index()
    )
    summary.to_csv(RESULTS_TABLES / "gc_distribution_summary_v0_2.csv", index=False)

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.violinplot(data=plot_df, x="candidate_type", y="gc_content", inner="box", cut=0, linewidth=1, ax=ax, color="#d9e6f2")
    ax.set_xlabel("candidate type")
    ax.set_ylabel("GC content")
    ax.set_title("v0.2 GC distribution by candidate type")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig2_gc_distribution_v0_2.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig3_single_landscape(pairs: pd.DataFrame, singles: pd.DataFrame) -> None:
    pair = resolve_pair(pairs, CONFIG["mutation_landscape_pair_id"])
    target = normalize_sequence(pair["target_dna"])
    df = singles.loc[singles["pair_id"] == pair["pair_id"]].copy()
    df["position"] = df["mutation_positions"].astype(int)
    df["alt_base"] = df["mutated_bases"]
    df["sequence_proxy_score"] = df.apply(
        lambda row: compute_sequence_metrics(target, row["candidate_dna"], tuple(CONFIG["sequence_baseline_k_values"]))["proxy_score"],
        axis=1,
    )
    summary = df.groupby(["position", "alt_base"], as_index=False)["sequence_proxy_score"].mean()
    summary.to_csv(RESULTS_TABLES / "single_mutation_sequence_proxy_landscape_v0_2.csv", index=False)

    positions = sorted(summary["position"].unique())
    matrix = pd.DataFrame(index=list("ACGT"), columns=positions, dtype=float)
    for pos in positions:
        for base in "ACGT":
            match = summary[(summary["position"] == pos) & (summary["alt_base"] == base)]
            matrix.loc[base, pos] = float(match["sequence_proxy_score"].iloc[0]) if not match.empty else np.nan

    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), height_ratios=[1, 1.4])
    mean_by_pos = summary.groupby("position")["sequence_proxy_score"].mean().reindex(positions)
    axes[0].plot(positions, mean_by_pos.values, marker="o", color="#dd8452", linewidth=1.8)
    axes[0].set_xlabel("single-mutation position")
    axes[0].set_ylabel("mean sequence-only proxy score")
    axes[0].set_title(f"Single-mutation sequence-proxy landscape: {pair['pdb_id']}")
    axes[0].set_ylim(0, 1)
    sns.heatmap(matrix, ax=axes[1], cmap="mako", vmin=0, vmax=1, cbar_kws={"label": "sequence-only proxy score"})
    axes[1].set_xlabel("single-mutation position")
    axes[1].set_ylabel("substituted base")
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig3_single_mutation_sequence_proxy_landscape_v0_2.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig4_double_landscape(pairs: pd.DataFrame, doubles: pd.DataFrame) -> None:
    pair = resolve_pair(pairs, CONFIG["mutation_landscape_pair_id"])
    target = normalize_sequence(pair["target_dna"])
    df = doubles.loc[doubles["pair_id"] == pair["pair_id"]].copy()
    positions = df["mutation_positions"].astype(str).str.split(";", n=1, expand=True)
    df["i"] = positions.iloc[:, 0].astype(int)
    df["j"] = positions.iloc[:, 1].astype(int)
    df["sequence_proxy_score"] = df.apply(
        lambda row: compute_sequence_metrics(target, row["candidate_dna"], tuple(CONFIG["sequence_baseline_k_values"]))["proxy_score"],
        axis=1,
    )
    summary = df.groupby(["i", "j"], as_index=False)["sequence_proxy_score"].mean()
    summary.to_csv(RESULTS_TABLES / "double_mutation_sequence_proxy_landscape_v0_2.csv", index=False)

    length = len(target)
    matrix = pd.DataFrame(np.nan, index=range(1, length + 1), columns=range(1, length + 1))
    for _, row in summary.iterrows():
        matrix.loc[int(row["i"]), int(row["j"])] = row["sequence_proxy_score"]

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(8, 7))
    mask = np.tril(np.ones_like(matrix, dtype=bool))
    sns.heatmap(
        matrix,
        mask=mask,
        ax=ax,
        cmap="viridis",
        vmin=0,
        vmax=1,
        cbar_kws={"label": "mean sequence-only proxy score"},
    )
    ax.set_xlabel("mutation position j")
    ax.set_ylabel("mutation position i")
    ax.set_title(f"Double-mutation sequence-proxy landscape: {pair['pdb_id']}")
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig4_double_mutation_sequence_proxy_landscape_v0_2.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig5_proxy_distribution(scored: pd.DataFrame) -> None:
    order = ["target", "single_mut", "double_mut", "gc_matched_random", "random_dna"]
    plot_df = scored[scored["candidate_type"].isin(order)].copy()
    plot_df["candidate_type"] = pd.Categorical(plot_df["candidate_type"], categories=order, ordered=True)
    summary = (
        plot_df.groupby("candidate_type", observed=False)["sequence_proxy_score"]
        .agg(["count", "mean", "std"])
        .reset_index()
    )
    summary.to_csv(RESULTS_TABLES / "sequence_proxy_distribution_summary_v0_2.csv", index=False)

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.violinplot(data=plot_df, x="candidate_type", y="sequence_proxy_score", inner=None, cut=0, linewidth=1, ax=ax, color="#cfd8dc")
    sns.boxplot(data=plot_df, x="candidate_type", y="sequence_proxy_score", width=0.18, ax=ax, showfliers=False, color="#4c72b0")
    ax.set_xlabel("candidate type")
    ax.set_ylabel("sequence-only proxy score")
    ax.set_title("v0.2 target and control sequence-only proxy distribution")
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig5_sequence_proxy_distribution_v0_2.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_figure_notes() -> None:
    notes = """# Figure Notes for Benchmark v0.2

Retrieval/audit date: 2026-08-29

## fig1_dataset_overview_v0_2.png

Can show: the size of the curated v0.2 perturbation benchmark per PDB-derived target.
Cannot show: binding specificity or off-target risk.
Category: dataset result / pipeline validation.

## fig2_gc_distribution_v0_2.png

Can show: GC-matched controls track target GC content more closely than fully random DNA.
Cannot show: binding preference or protein-conditioned specificity.
Category: dataset result / pipeline validation.

## fig3_single_mutation_sequence_proxy_landscape_v0_2.png

Can show: how sequence-only proxy metrics change under single substitutions for one curated target.
Cannot show: biological specificity landscape or mutation effect on binding.
Category: sequence-only baseline.

## fig4_double_mutation_sequence_proxy_landscape_v0_2.png

Can show: how sequence-only proxy metrics change under double substitutions for one curated target.
Cannot show: epistasis, binding energy, or off-target risk.
Category: sequence-only baseline.

## fig5_sequence_proxy_distribution_v0_2.png

Can show: separation induced by sequence similarity between target, mutants, GC controls, and random controls.
Cannot show: protein-DNA binding specificity or calibrated safety margin.
Category: sequence-only baseline / pipeline sanity check.

## Shared Limitation

No v0.2 figure should be interpreted as a biological specificity result. The current scoring is not protein-conditioned and has no quantitative experimental specificity ground truth attached.
"""
    (RESULTS_FIGURES / "FIGURE_NOTES.md").write_text(notes, encoding="utf-8")


def main() -> None:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs_v0_2.csv")
    benchmark = pd.read_csv(PROCESSED_DIR / "benchmark_v0_2.csv", low_memory=False)
    scored = pd.read_csv(PROCESSED_DIR / "benchmark_v0_2_scored.csv", low_memory=False)
    singles = pd.read_csv(PROCESSED_DIR / "single_mutants_v0_2.csv")
    doubles = pd.read_csv(PROCESSED_DIR / "double_mutants_v0_2.csv")
    fig1_dataset_overview(benchmark, pairs)
    fig2_gc_distribution(scored)
    fig3_single_landscape(pairs, singles)
    fig4_double_landscape(pairs, doubles)
    fig5_proxy_distribution(scored)
    write_figure_notes()
    print(f"Saved v0.2 figures to {RESULTS_FIGURES}")


if __name__ == "__main__":
    main()
