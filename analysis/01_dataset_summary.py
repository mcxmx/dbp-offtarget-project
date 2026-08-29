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

from src.utils import ensure_dir, load_yaml, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_FIGURES = ensure_dir(ROOT / "results" / "figures")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")


def main() -> None:
    benchmark = pd.read_csv(PROCESSED_DIR / "benchmark_v0_1.csv")
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")

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
    counts.to_csv(RESULTS_TABLES / "dataset_overview_counts.csv")

    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [3.2, 1.2]})

    log_counts = np.log10(counts + 1)
    sns.heatmap(
        log_counts,
        ax=axes[0],
        cmap="viridis",
        cbar_kws={"label": "log10(count + 1)"},
        linewidths=0.3,
        linecolor="white",
    )
    axes[0].set_xlabel("candidate type")
    axes[0].set_ylabel("pair_id")
    axes[0].set_title("Dataset overview by pair")

    lengths = pairs.set_index("pair_id").loc[order, "dna_length"]
    axes[1].barh(range(len(order)), lengths.values, color="#4c72b0")
    axes[1].set_yticks(range(len(order)))
    axes[1].set_yticklabels([])
    axes[1].set_xlabel("DNA length")
    axes[1].set_title("Target length")
    axes[1].invert_yaxis()

    fig.tight_layout()
    fig.savefig(RESULTS_FIGURES / "fig1_dataset_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {RESULTS_FIGURES / 'fig1_dataset_overview.png'}")


if __name__ == "__main__":
    main()
