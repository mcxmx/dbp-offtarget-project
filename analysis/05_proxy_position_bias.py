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


def kmers_covering_position(length: int, position: int, k: int) -> int:
    start_min = max(1, position - k + 1)
    start_max = min(position, length - k + 1)
    return max(0, start_max - start_min + 1)


def build_summary(pair: pd.Series, singles: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    target = normalize_sequence(pair["target_dna"])
    pair_id = pair["pair_id"]
    df = singles.loc[singles["pair_id"] == pair_id].copy()
    if df.empty:
        raise ValueError(f"No single mutants found for {pair_id}")
    df["position"] = df["mutation_positions"].astype(int)
    rows = []
    for _, row in df.iterrows():
        metrics = compute_sequence_metrics(target, row["candidate_dna"], tuple(CONFIG["sequence_baseline_k_values"]))
        length = len(target)
        rows.append(
            {
                "pair_id": pair_id,
                "pdb_id": pair["pdb_id"],
                "position": int(row["position"]),
                "mutated_bases": row["mutated_bases"],
                "hamming_similarity": 1.0 - metrics["hamming_distance"] / length,
                "edit_similarity": 1.0 - metrics["edit_distance"] / length,
                "kmer3_jaccard": metrics["kmer3_jaccard"],
                "kmer4_jaccard": metrics["kmer4_jaccard"],
                "combined_sequence_proxy": metrics["proxy_score"],
            }
        )
    long = pd.DataFrame(rows).melt(
        id_vars=["pair_id", "pdb_id", "position", "mutated_bases"],
        value_vars=[
            "hamming_similarity",
            "edit_similarity",
            "kmer3_jaccard",
            "kmer4_jaccard",
            "combined_sequence_proxy",
        ],
        var_name="metric",
        value_name="sequence_only_proxy_value",
    )
    summary = long.groupby(["pair_id", "pdb_id", "position", "metric"], as_index=False)["sequence_only_proxy_value"].mean()
    artifact_rows = []
    for position in range(1, len(target) + 1):
        artifact_rows.append(
            {
                "pair_id": pair_id,
                "pdb_id": pair["pdb_id"],
                "position": position,
                "length": len(target),
                "kmer3_windows_covering_position": kmers_covering_position(len(target), position, 3),
                "kmer4_windows_covering_position": kmers_covering_position(len(target), position, 4),
            }
        )
    artifact = pd.DataFrame(artifact_rows)
    return summary, artifact


def plot(summary: pd.DataFrame, artifact: pd.DataFrame, pair_id: str) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, height_ratios=[2.2, 1])
    order = [
        "hamming_similarity",
        "edit_similarity",
        "kmer3_jaccard",
        "kmer4_jaccard",
        "combined_sequence_proxy",
    ]
    plot_df = summary.copy()
    plot_df["metric"] = pd.Categorical(plot_df["metric"], categories=order, ordered=True)
    sns.lineplot(
        data=plot_df,
        x="position",
        y="sequence_only_proxy_value",
        hue="metric",
        marker="o",
        ax=axes[0],
    )
    axes[0].set_ylabel("mean sequence-only proxy value")
    axes[0].set_title(f"Position effect of sequence-only proxy metrics: {pair_id}")
    axes[0].set_ylim(0, 1.02)
    axes[0].legend(title="metric", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)

    axes[1].plot(
        artifact["position"],
        artifact["kmer3_windows_covering_position"],
        marker="o",
        label="3-mer windows covering position",
        color="#4c72b0",
    )
    axes[1].plot(
        artifact["position"],
        artifact["kmer4_windows_covering_position"],
        marker="s",
        label="4-mer windows covering position",
        color="#dd8452",
    )
    axes[1].set_xlabel("single-mutation position")
    axes[1].set_ylabel("overlapping k-mers")
    axes[1].legend(title=None)
    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig_proxy_position_bias.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs_v0_2.csv")
    singles = pd.read_csv(PROCESSED_DIR / "single_mutants_v0_2.csv")
    pair = resolve_pair(pairs, CONFIG["mutation_landscape_pair_id"])
    summary, artifact = build_summary(pair, singles)
    summary.to_csv(RESULTS_TABLES / "proxy_position_bias_summary.csv", index=False)
    artifact.to_csv(RESULTS_TABLES / "proxy_position_bias_kmer_windows.csv", index=False)
    metric_ranges = (
        summary.groupby("metric")["sequence_only_proxy_value"]
        .agg(["min", "max", "mean"])
        .reset_index()
    )
    metric_ranges["range"] = metric_ranges["max"] - metric_ranges["min"]
    metric_ranges.to_csv(RESULTS_TABLES / "proxy_position_bias_metric_ranges.csv", index=False)
    plot(summary, artifact, pair["pair_id"])
    print(f"Saved {RESULTS_FIGURES / 'fig_proxy_position_bias.png'}")
    print(metric_ranges.to_string(index=False))


if __name__ == "__main__":
    main()

