from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
RESULTS = ensure_dir(ROOT / "results" / "v0_4")
TABLES = ensure_dir(RESULTS / "tables")
FIGURES = ensure_dir(RESULTS / "figures")

PROTEIN_ORDER = ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]
BASELINE_LABELS = {
    "sequence_hamming": "Seq Hamming",
    "sequence_edit": "Seq edit",
    "sequence_kmer3": "Seq 3-mer",
    "sequence_kmer4": "Seq 4-mer",
    "NA-MPNN_structural_ppm": "NA-MPNN PPM",
    "DeepPBS": "DeepPBS",
    "SimpleProteinConditionalBaseline": "Simple PC",
}


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def savefig(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def fig1_macro_overview(macro: pd.DataFrame, noise: pd.DataFrame) -> None:
    sub = macro[macro["metric"] == "spearman"].copy()
    order = [
        "sequence_hamming",
        "sequence_edit",
        "sequence_kmer3",
        "sequence_kmer4",
        "SimpleProteinConditionalBaseline",
        "DeepPBS",
        "NA-MPNN_structural_ppm",
    ]
    sub["baseline"] = pd.Categorical(sub["baseline"], categories=order, ordered=True)
    sub = sub.sort_values("baseline")
    replicate_ref = float(noise.loc[noise["score_type"] == "e_score", "spearman_correlation"].median())

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    colors = ["#4c78a8" if n > 0 else "#c7c7c7" for n in sub["n_proteins"]]
    heights = sub["median"].fillna(0.0).to_numpy()
    bars = ax.bar([BASELINE_LABELS.get(x, x) for x in sub["baseline"]], heights, color=colors, edgecolor="#444444")
    for bar, median, n in zip(bars, sub["median"], sub["n_proteins"]):
        label = "NA" if pd.isna(median) else f"{median:.2f}"
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02, f"{label}\nN={int(n)}", ha="center", va="bottom", fontsize=8)
        if int(n) == 0:
            bar.set_hatch("//")
    ax.axhline(replicate_ref, color="#d62728", linestyle="--", linewidth=1.4, label=f"Replicate ref. {replicate_ref:.2f}")
    ax.set_ylabel("Macro median Spearman")
    ax.set_ylim(-0.15, max(0.75, np.nanmax(sub["median"].to_numpy(dtype=float)) + 0.15))
    ax.set_title("v0.4 baseline ranking performance overview")
    ax.legend(frameon=False, loc="upper right")
    ax.tick_params(axis="x", rotation=25)
    savefig(FIGURES / "fig_v0_4_1_baseline_performance_overview.png")


def fig2_per_protein_heatmap(per_protein: pd.DataFrame) -> None:
    baselines = [
        "sequence_hamming",
        "sequence_edit",
        "sequence_kmer3",
        "sequence_kmer4",
        "NA-MPNN_structural_ppm",
        "DeepPBS",
        "SimpleProteinConditionalBaseline",
    ]
    pivot = per_protein.pivot(index="baseline", columns="protein_id", values="spearman").reindex(baselines)[PROTEIN_ORDER]
    labels = pivot.apply(lambda col: col.map(lambda x: "NA" if pd.isna(x) else f"{x:.2f}"))
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="vlag",
        center=0.0,
        vmin=-0.4,
        vmax=0.4,
        annot=labels,
        fmt="",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Spearman"},
    )
    for i, baseline in enumerate(pivot.index):
        for j, protein_id in enumerate(pivot.columns):
            if pd.isna(pivot.loc[baseline, protein_id]):
                ax.text(j + 0.5, i + 0.5, "NA", ha="center", va="center", color="#666666", fontsize=9)
    ax.set_yticklabels([BASELINE_LABELS.get(x.get_text(), x.get_text()) for x in ax.get_yticklabels()], rotation=0)
    ax.set_xlabel("Designed DBP")
    ax.set_ylabel("Baseline")
    ax.set_title("Per-protein Spearman on PBM E-score ranking")
    savefig(FIGURES / "fig_v0_4_2_per_protein_spearman_heatmap.png")


def _spearman_for(df: pd.DataFrame, score_col: str) -> float:
    sub = df[["experimental_escore_consensus", score_col]].dropna()
    if sub.empty or sub[score_col].nunique() < 2:
        return np.nan
    return float(sub["experimental_escore_consensus"].corr(sub[score_col], method="spearman"))


def fig3_prediction_examples(scored: pd.DataFrame) -> None:
    proteins = ["DBP35", "DBP48"]
    score_cols = [("sequence_kmer3_score", "Seq 3-mer"), ("structural_ppm_score", "NA-MPNN PPM")]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.5), sharex=False, sharey=False)
    for row_idx, protein_id in enumerate(proteins):
        g = scored[scored["protein_id"] == protein_id]
        for col_idx, (score_col, label) in enumerate(score_cols):
            ax = axes[row_idx, col_idx]
            sub = g[["experimental_escore_consensus", score_col]].dropna()
            if sub.empty:
                ax.text(0.5, 0.5, "Not evaluable", transform=ax.transAxes, ha="center", va="center")
                ax.set_axis_off()
                continue
            sns.scatterplot(
                data=sub,
                x=score_col,
                y="experimental_escore_consensus",
                s=7,
                alpha=0.35,
                linewidth=0,
                color="#2f6f8f",
                ax=ax,
            )
            rho = _spearman_for(g, score_col)
            ax.set_title(f"{protein_id}: {label}, Spearman {rho:.2f}")
            ax.set_xlabel("Baseline score")
            ax.set_ylabel("Processed PBM E-score consensus")
    savefig(FIGURES / "fig_v0_4_3_prediction_vs_experimental_examples.png")


def fig4_vs_replicate_reference(per_protein: pd.DataFrame, noise: pd.DataFrame) -> None:
    use = per_protein[per_protein["baseline"].isin(["sequence_kmer3", "NA-MPNN_structural_ppm"])].copy()
    use["baseline_label"] = use["baseline"].map(BASELINE_LABELS)
    use["protein_id"] = pd.Categorical(use["protein_id"], categories=PROTEIN_ORDER, ordered=True)
    ref = noise[noise["score_type"] == "e_score"][["protein_id", "spearman_correlation"]].copy()
    ref["protein_id"] = pd.Categorical(ref["protein_id"], categories=PROTEIN_ORDER, ordered=True)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    sns.barplot(data=use, x="protein_id", y="spearman", hue="baseline_label", ax=ax, palette=["#4c78a8", "#f58518"])
    sns.scatterplot(data=ref, x="protein_id", y="spearman_correlation", ax=ax, color="#d62728", s=55, label="Replicate ref.")
    ax.axhline(ref["spearman_correlation"].median(), color="#d62728", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.set_ylabel("Spearman")
    ax.set_xlabel("Designed DBP")
    ax.set_ylim(-0.45, 0.85)
    ax.set_title("Baseline performance relative to empirical replicate agreement")
    ax.legend(frameon=False, loc="upper left")
    savefig(FIGURES / "fig_v0_4_4_performance_vs_replicate_reference.png")


def fig5_by_distance(distance: pd.DataFrame) -> None:
    use = distance[
        (distance["baseline"].isin(["sequence_kmer3", "sequence_kmer4", "NA-MPNN_structural_ppm"]))
        & (distance["evaluation_status"] == "evaluated")
    ].copy()
    summary = (
        use.groupby(["baseline", "motif_distance_bin"], as_index=False)
        .agg(median_spearman=("spearman", "median"), n_groups=("spearman", "count"))
    )
    distance_order = ["0", "1", "2", "3+"]
    summary["motif_distance_bin"] = pd.Categorical(summary["motif_distance_bin"], categories=distance_order, ordered=True)
    summary = summary.sort_values("motif_distance_bin")
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    sns.lineplot(
        data=summary,
        x="motif_distance_bin",
        y="median_spearman",
        hue="baseline",
        marker="o",
        ax=ax,
        palette=["#4c78a8", "#72b7b2", "#f58518"],
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, [BASELINE_LABELS.get(label, label) for label in labels], frameon=False, title="")
    ax.set_xlabel("RC-aware motif Hamming distance bin")
    ax.set_ylabel("Median Spearman within bin")
    ax.set_title("Performance by motif-distance regime")
    ax.set_ylim(-0.45, 0.65)
    savefig(FIGURES / "fig_v0_4_5_performance_by_motif_distance.png")


def fig6_failure_landscape(failures: pd.DataFrame, resolution: pd.DataFrame) -> None:
    counts = failures.groupby(["protein_id", "failure_category"]).size().reset_index(name="n_cases")
    counts["failure_label"] = counts["failure_category"].map(
        {
            "experimental_high_sequence_proxy_low": "Experimental high / sequence low",
            "sequence_proxy_high_experimental_low": "Sequence high / experimental low",
        }
    )
    counts["protein_id"] = pd.Categorical(counts["protein_id"], categories=PROTEIN_ORDER, ordered=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    sns.barplot(data=counts, x="protein_id", y="n_cases", hue="failure_label", ax=axes[0], palette=["#e45756", "#72b7b2"])
    axes[0].set_xlabel("Designed DBP")
    axes[0].set_ylabel("Sequence-vs-experiment cases")
    axes[0].set_title("Failure-case categories")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].legend(frameon=False, title="")

    res = resolution[resolution["baseline"] == "NA-MPNN_structural_ppm"].copy()
    res["protein_id"] = pd.Categorical(res["protein_id"], categories=PROTEIN_ORDER, ordered=True)
    sns.barplot(data=res, x="protein_id", y="n_v0_3_1_disagreement_candidates", color="#d0d0d0", ax=axes[1], label="Total / not-evaluable background")
    sns.barplot(data=res, x="protein_id", y="n_resolved", color="#4c78a8", ax=axes[1], label="Resolved by NA-MPNN top 10%")
    axes[1].set_xlabel("Designed DBP")
    axes[1].set_ylabel("v0.3.1 disagreement candidates")
    axes[1].set_title("Resolution of sequence-vs-experiment disagreements")
    axes[1].tick_params(axis="x", rotation=25)
    axes[1].legend(frameon=False)
    savefig(FIGURES / "fig_v0_4_6_failure_case_landscape.png")


def write_figure_notes() -> None:
    notes = """# v0.4 Figure Notes

All v0.4 panels evaluate processed experimental uPBM E-score ranking, not binding affinity or off-target prediction.

- `fig_v0_4_1_baseline_performance_overview.png`: macro median per-protein Spearman. This shows broad ranking performance and evaluation coverage; it does not prove absolute calibration.
- `fig_v0_4_2_per_protein_spearman_heatmap.png`: per-DBP Spearman. Gray/NA cells are not-evaluable methods, not zero-valued failures.
- `fig_v0_4_3_prediction_vs_experimental_examples.png`: DBP35 and DBP48 examples for sequence-only and NA-MPNN structural PPM scores. These are diagnostic mappings to PBM 7-mer ranking.
- `fig_v0_4_4_performance_vs_replicate_reference.png`: compares evaluated baselines with empirical replicate agreement. Replicate agreement is an assay reproducibility reference, not a strict upper bound.
- `fig_v0_4_5_performance_by_motif_distance.png`: checks whether performance changes with RC-aware motif-distance regime.
- `fig_v0_4_6_failure_case_landscape.png`: separates total sequence-vs-experiment disagreement cases from examples resolved by NA-MPNN where predictions exist.
"""
    (FIGURES / "FIGURE_NOTES.md").write_text(notes, encoding="utf-8")


def main() -> None:
    setup_style()
    macro = pd.read_csv(TABLES / "baseline_performance_macro.csv")
    per_protein = pd.read_csv(TABLES / "baseline_performance_per_protein.csv")
    scored = pd.read_parquet(ROOT / "data" / "processed" / "v0_4" / "v0_4_scored_candidates.parquet")
    distance = pd.read_csv(TABLES / "performance_by_sequence_distance.csv")
    failures = pd.read_parquet(TABLES / "baseline_failure_cases.parquet")
    resolution = pd.read_csv(TABLES / "failure_resolution_summary.csv")
    noise = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "experimental_noise_ceiling.csv")

    fig1_macro_overview(macro, noise)
    fig2_per_protein_heatmap(per_protein)
    fig3_prediction_examples(scored)
    fig4_vs_replicate_reference(per_protein, noise)
    fig5_by_distance(distance)
    fig6_failure_landscape(failures, resolution)
    write_figure_notes()
    print(f"Wrote v0.4 figures to {FIGURES.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
