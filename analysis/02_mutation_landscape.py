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
RESULTS_FIGURES = ensure_dir(ROOT / "results" / "figures")
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


def single_landscape(pair_id: str, pairs: pd.DataFrame, singles: pd.DataFrame) -> pd.DataFrame:
    pair = resolve_pair(pairs, pair_id)
    resolved_pair_id = pair["pair_id"]
    target = normalize_sequence(pair["target_dna"])
    df = singles.loc[singles["pair_id"] == resolved_pair_id].copy()
    if df.empty:
        raise ValueError(f"No single mutants found for {resolved_pair_id}")
    df["position"] = df["mutation_positions"].astype(str).str.split(";", n=1, expand=True).iloc[:, 0].astype(int)
    df["alt_base"] = df["mutated_bases"]
    df["proxy_score"] = df.apply(lambda r: compute_sequence_metrics(target, r["candidate_dna"], tuple(CONFIG["sequence_baseline_k_values"]))["proxy_score"], axis=1)
    summary = df.groupby(["position", "alt_base"])["proxy_score"].mean().reset_index()
    summary.to_csv(RESULTS_TABLES / "single_mutation_landscape_summary.csv", index=False)
    return summary, target


def double_landscape(pair_id: str, pairs: pd.DataFrame, doubles: pd.DataFrame) -> pd.DataFrame:
    pair = resolve_pair(pairs, pair_id)
    resolved_pair_id = pair["pair_id"]
    target = normalize_sequence(pair["target_dna"])
    df = doubles.loc[doubles["pair_id"] == resolved_pair_id].copy()
    if df.empty:
        raise ValueError(f"No double mutants found for {resolved_pair_id}")
    positions = df["mutation_positions"].astype(str).str.split(";", n=1, expand=True)
    if positions.shape[1] != 2:
        raise ValueError(f"Unexpected mutation_positions format for {resolved_pair_id}")
    df["i"] = positions.iloc[:, 0].astype(int)
    df["j"] = positions.iloc[:, 1].astype(int)
    df["proxy_score"] = df.apply(lambda r: compute_sequence_metrics(target, r["candidate_dna"], tuple(CONFIG["sequence_baseline_k_values"]))["proxy_score"], axis=1)
    summary = df.groupby(["i", "j"])["proxy_score"].mean().reset_index()
    summary.to_csv(RESULTS_TABLES / "double_mutation_landscape_summary.csv", index=False)
    return summary, target


def plot_single(summary: pd.DataFrame, target: str, pair_id: str) -> None:
    positions = sorted(summary["position"].unique())
    bases = list("ACGT")
    matrix = pd.DataFrame(index=bases, columns=positions, dtype=float)
    for pos in positions:
        for base in bases:
            match = summary[(summary["position"] == pos) & (summary["alt_base"] == base)]
            matrix.loc[base, pos] = float(match["proxy_score"].iloc[0]) if not match.empty else np.nan
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), height_ratios=[1, 1.3])
    mean_by_pos = summary.groupby("position")["proxy_score"].mean().reindex(positions)
    axes[0].plot(positions, mean_by_pos.values, marker="o", color="#dd8452", linewidth=1.8)
    axes[0].set_xlabel("mutation position")
    axes[0].set_ylabel("mean proxy score")
    axes[0].set_title(f"Single mutation landscape: {pair_id}")
    axes[0].set_ylim(0, 1)
    sns.heatmap(matrix, ax=axes[1], cmap="mako", vmin=0, vmax=1, cbar_kws={"label": "proxy score"})
    axes[1].set_xlabel("mutation position")
    axes[1].set_ylabel("substituted base")
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig3_single_mutation_landscape.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_double(summary: pd.DataFrame, pair_id: str, target: str) -> None:
    L = len(target)
    matrix = pd.DataFrame(np.nan, index=range(1, L + 1), columns=range(1, L + 1))
    for _, row in summary.iterrows():
        i = int(row["i"])
        j = int(row["j"])
        matrix.loc[i, j] = row["proxy_score"]
    mask = np.tril(np.ones_like(matrix, dtype=bool))
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(9, 8))
    sns.heatmap(
        matrix,
        mask=mask,
        ax=ax,
        cmap="viridis",
        vmin=0,
        vmax=1,
        cbar_kws={"label": "mean proxy score"},
    )
    ax.set_xlabel("mutation position j")
    ax.set_ylabel("mutation position i")
    ax.set_title(f"Double mutation landscape: {pair_id}")
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig4_double_mutation_landscape.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")
    singles = pd.read_csv(PROCESSED_DIR / "single_mutants.csv")
    doubles = pd.read_csv(PROCESSED_DIR / "double_mutants.csv")
    pair_id = CONFIG["mutation_landscape_pair_id"]
    single_summary, target = single_landscape(pair_id, pairs, singles)
    double_summary, _ = double_landscape(pair_id, pairs, doubles)
    plot_single(single_summary, target, pair_id)
    plot_double(double_summary, pair_id, target)
    print(f"Saved fig3 and fig4 for {pair_id}")


if __name__ == "__main__":
    main()
