from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.utils import ensure_dir, project_root, sequence_identity
from src.v0_4_evaluation import compute_ranking_metrics


ROOT = project_root()
PROCESSED_V031 = ROOT / "data" / "processed" / "v0_3_1"
V042_RESULTS = ensure_dir(ROOT / "results" / "v0_4_2")
V042_TABLES = ensure_dir(V042_RESULTS / "tables")
V042_FIGURES = ensure_dir(V042_RESULTS / "figures")
V042_DATA = ensure_dir(ROOT / "data" / "processed" / "v0_4_2")
V042_DOCS = ensure_dir(ROOT / "docs" / "v0_4_2")
SEED = 42
DESIGNED_PROTEINS = ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]


def load_designed() -> pd.DataFrame:
    df = pd.read_parquet(PROCESSED_V031 / "designed_dbp_upbm_rc_class_v0_3_1.parquet").copy()
    df = df.rename(columns={"experimental_escore_consensus": "experimental_score"})
    df["experimental_percentile"] = df.groupby("protein_id")["experimental_score"].rank(
        pct=True, ascending=True
    )
    return df


def load_predictions(designed: pd.DataFrame) -> pd.DataFrame:
    seq = pd.read_parquet(
        PROCESSED_V031 / "designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1.parquet"
    )[
        [
            "protein_id",
            "canonical_7mer",
            "kmer3_jaccard_to_paper_motif_rc_aware",
            "hamming_similarity_to_paper_motif_rc_aware",
        ]
    ].rename(
        columns={
            "kmer3_jaccard_to_paper_motif_rc_aware": "sequence_score",
            "hamming_similarity_to_paper_motif_rc_aware": "sequence_hamming_score",
        }
    )
    simple = pd.read_parquet(
        ROOT / "results" / "v0_4_1" / "tables" / "simple_pc_designed_predictions.parquet"
    )[["protein_id", "canonical_rc", "simple_pc_score"]].rename(
        columns={"canonical_rc": "canonical_7mer"}
    )
    frozen = pd.read_parquet(
        ROOT / "results" / "v0_4_2" / "tables" / "frozen_plm_designed_predictions.parquet"
    )[["protein_id", "canonical_rc", "frozen_plm_score"]].rename(
        columns={"canonical_rc": "canonical_7mer"}
    )
    nampnn_path = ROOT / "results" / "v0_4" / "tables" / "nampnn_predictions.parquet"
    if nampnn_path.exists():
        nampnn = pd.read_parquet(nampnn_path)[
            ["protein_id", "canonical_7mer", "prediction_score"]
        ].rename(columns={"prediction_score": "nampnn_score"})
    else:
        nampnn = pd.DataFrame(columns=["protein_id", "canonical_7mer", "nampnn_score"])

    merged = designed.merge(seq, on=["protein_id", "canonical_7mer"], how="left")
    merged = merged.merge(simple, on=["protein_id", "canonical_7mer"], how="left")
    merged = merged.merge(frozen, on=["protein_id", "canonical_7mer"], how="left")
    merged = merged.merge(nampnn, on=["protein_id", "canonical_7mer"], how="left")
    merged["deeppbs_score"] = np.nan
    for score_col in [
        "sequence_score",
        "simple_pc_score",
        "frozen_plm_score",
        "nampnn_score",
        "deeppbs_score",
    ]:
        pct_col = f"{score_col}_percentile"
        merged[pct_col] = merged.groupby("protein_id")[score_col].rank(
            pct=True, ascending=True
        )
    return merged


def best_motif_distance(candidate: str, motif: str) -> int:
    candidate = str(candidate).upper()
    motif = str(motif).upper()
    if not candidate or not motif or len(candidate) < len(motif):
        return np.nan
    candidate_orientations = {candidate, reverse_complement(candidate)}
    motif_orientations = {motif, reverse_complement(motif)}
    distances = []
    for oriented_candidate in candidate_orientations:
        for start in range(len(oriented_candidate) - len(motif) + 1):
            window = oriented_candidate[start : start + len(motif)]
            for oriented_motif in motif_orientations:
                distances.append(sum(a != b for a, b in zip(window, oriented_motif)))
    return int(min(distances))


def add_descriptive_features(df: pd.DataFrame) -> pd.DataFrame:
    motifs = (
        df[["protein_id", "designed_binding_site_motif"]]
        .drop_duplicates("protein_id")
        .set_index("protein_id")["designed_binding_site_motif"]
        .to_dict()
    )
    df = df.copy()
    df["motif_hamming_distance"] = [
        best_motif_distance(sequence, motifs[protein])
        for protein, sequence in zip(df["protein_id"], df["canonical_7mer"])
    ]
    df["motif_distance_bin"] = df["motif_hamming_distance"].map(
        lambda value: str(int(value)) if pd.notna(value) and value < 3 else "3+"
    )
    return df


def disagreement_mask(df: pd.DataFrame) -> pd.Series:
    counts = pd.read_csv(
        ROOT / "results" / "v0_3_1" / "tables" / "all_disagreement_candidate_counts.csv"
    ).set_index("protein_id")
    mask = pd.Series(False, index=df.index)
    for protein_id, group in df.groupby("protein_id"):
        thresholds = counts.loc[protein_id]
        local = (
            (group["experimental_score"] >= float(thresholds["experimental_score_threshold"]))
            & (
                group["sequence_hamming_score"]
                <= float(thresholds["sequence_similarity_threshold"]) + 1e-12
            )
        )
        mask.loc[group.index] = local
    return mask


def write_disagreement_resolution(df: pd.DataFrame) -> pd.DataFrame:
    candidate_mask = disagreement_mask(df)
    method_cols = {
        "sequence_kmer3": "sequence_score",
        "SimpleProteinConditionalBaseline": "simple_pc_score",
        "FrozenPLMProteinConditionalBaseline": "frozen_plm_score",
        "NA-MPNN diagnostic": "nampnn_score",
        "DeepPBS": "deeppbs_score",
    }
    rows = []
    for protein_id, group in df.groupby("protein_id", sort=True):
        candidates = group.loc[candidate_mask.loc[group.index]]
        for method, score_col in method_cols.items():
            evaluable = candidates[score_col].notna()
            resolved = evaluable & (
                candidates[f"{score_col}_percentile"] >= 0.90
            )
            rows.append(
                {
                    "protein_id": protein_id,
                    "method": method,
                    "n_total_candidates": int(len(candidates)),
                    "n_evaluable": int(evaluable.sum()),
                    "n_resolved": int(resolved.sum()),
                    "n_unresolved": int((evaluable & ~resolved).sum()),
                    "resolution_rate_among_evaluable": (
                        float(resolved.sum() / evaluable.sum())
                        if evaluable.sum()
                        else np.nan
                    ),
                    "resolution_rate_among_all_candidates": float(
                        resolved.sum() / len(candidates)
                    )
                    if len(candidates)
                    else np.nan,
                    "threshold_definition": (
                        "v0.3.1 candidate set; method prediction percentile >= 0.90 "
                        "within protein"
                    ),
                    "evaluation_status": (
                        "evaluated"
                        if evaluable.any()
                        else "not_evaluable_missing_prediction"
                    ),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(
        V042_TABLES / "disagreement_resolution_v0_4_2.csv", index=False
    )
    return out


def write_common_hard_set(df: pd.DataFrame) -> pd.DataFrame:
    core = ["sequence_score", "simple_pc_score", "frozen_plm_score"]
    high_experiment = df["experimental_percentile"] >= 0.95
    low_all_core = df[[f"{col}_percentile" for col in core]].le(0.50).all(axis=1)
    low_experiment = df["experimental_percentile"] <= 0.50
    high_all_core = df[[f"{col}_percentile" for col in core]].ge(0.95).all(axis=1)
    high_mask = high_experiment & low_all_core
    low_mask = low_experiment & high_all_core
    selected = df.loc[high_mask | low_mask].copy()
    selected["failure_type"] = np.where(
        high_mask.loc[selected.index],
        "high_experiment_low_all_core",
        "low_experiment_high_all_core",
    )
    selected["sequence_distance_to_target"] = selected["motif_hamming_distance"]
    selected["notes"] = (
        "Core common-hard diagnostic using sequence kmer3, SimplePC, and "
        "FrozenPLM; DeepPBS unavailable and NA-MPNN partial."
    )
    columns = [
        "protein_id",
        "canonical_7mer",
        "experimental_score",
        "experimental_percentile",
        "sequence_score",
        "sequence_score_percentile",
        "simple_pc_score",
        "simple_pc_score_percentile",
        "frozen_plm_score",
        "frozen_plm_score_percentile",
        "deeppbs_score",
        "nampnn_score",
        "failure_type",
        "sequence_distance_to_target",
        "motif_distance_bin",
        "notes",
    ]
    selected[columns].to_parquet(
        V042_DATA / "common_hard_specificity_cases.parquet", index=False
    )
    summary = (
        selected.groupby(["protein_id", "failure_type"], as_index=False)
        .size()
        .rename(columns={"size": "n_cases"})
    )
    summary.to_csv(V042_TABLES / "common_hard_case_counts.csv", index=False)
    return selected


def score_distribution_features(group: pd.DataFrame) -> dict[str, float]:
    values = group["experimental_score"].to_numpy(dtype=float)
    hist, _ = np.histogram(values, bins=10)
    probs = hist[hist > 0] / len(values)
    histogram_entropy = float(-(probs * np.log(probs)).sum() / np.log(10))
    shifted = values - np.min(values) + 1e-9
    top_n = max(1, int(math.ceil(len(values) * 0.01)))
    top_values = np.sort(values)[-top_n:]
    return {
        "experimental_score_std": float(np.std(values, ddof=1)),
        "experimental_score_variance": float(np.var(values, ddof=1)),
        "experimental_score_iqr": float(np.percentile(values, 75) - np.percentile(values, 25)),
        "experimental_score_q95": float(np.percentile(values, 95)),
        "experimental_score_q99": float(np.percentile(values, 99)),
        "top_1pct_mean": float(np.mean(top_values)),
        "top_1pct_minus_median": float(np.mean(top_values) - np.median(values)),
        "top_1pct_shifted_mass_fraction": float(np.sum(top_values - np.min(values) + 1e-9) / np.sum(shifted)),
        "score_histogram_entropy_normalized": histogram_entropy,
        "n_unique_experimental_scores": int(group["experimental_score"].nunique()),
    }


def write_performance_by_distance(df: pd.DataFrame) -> pd.DataFrame:
    methods = {
        "sequence_kmer3": "sequence_score",
        "SimpleProteinConditionalBaseline": "simple_pc_score",
        "FrozenPLMProteinConditionalBaseline": "frozen_plm_score",
        "NA-MPNN diagnostic": "nampnn_score",
        "DeepPBS": "deeppbs_score",
    }
    rows = []
    for (protein_id, distance_bin), group in df.groupby(
        ["protein_id", "motif_distance_bin"], sort=True
    ):
        for method, score_col in methods.items():
            sub = group[["experimental_score", score_col, "canonical_7mer"]].dropna()
            if len(sub) < 20 or sub[score_col].nunique() < 2:
                rows.append(
                    {
                        "protein_id": protein_id,
                        "motif_distance_bin": distance_bin,
                        "method": method,
                        "spearman": np.nan,
                        "ndcg_5pct": np.nan,
                        "n_rc_classes": int(len(sub)),
                        "evaluation_status": "not_evaluable_or_too_few_scores",
                    }
                )
                continue
            metrics = compute_ranking_metrics(
                sub, "experimental_score", score_col, seed=SEED
            )
            rows.append(
                {
                    "protein_id": protein_id,
                    "motif_distance_bin": distance_bin,
                    "method": method,
                    "spearman": metrics.spearman,
                    "ndcg_5pct": metrics.ndcg_5pct,
                    "n_rc_classes": metrics.n_rc_classes,
                    "evaluation_status": "evaluated",
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(V042_TABLES / "performance_by_motif_distance_v0_4_2.csv", index=False)
    return out


def write_difficulty_factors(df: pd.DataFrame) -> pd.DataFrame:
    groups = pd.read_csv(ROOT / "metadata" / "v0_3_1" / "designed_dbp_target_groups.csv")
    seq_perf = pd.read_csv(
        ROOT / "results" / "v0_4" / "tables" / "baseline_performance_per_protein.csv"
    )
    seq_perf = seq_perf.query("baseline == 'sequence_kmer3'")[
        ["protein_id", "spearman"]
    ].rename(columns={"spearman": "sequence_kmer3_spearman"})
    pc_perf = pd.read_csv(
        ROOT / "results" / "v0_4_1" / "tables" / "simple_pc_performance.csv"
    )
    pc_perf = pc_perf.query("dataset == 'designed_external'")[
        ["protein_id", "spearman"]
    ].rename(columns={"spearman": "simple_pc_spearman"})
    plm_perf = pd.read_csv(
        ROOT / "results" / "v0_4_2" / "tables" / "frozen_plm_performance.csv"
    )
    plm_perf = plm_perf.query("dataset == 'designed_external'")[
        ["protein_id", "spearman"]
    ].rename(columns={"spearman": "frozen_plm_spearman"})

    natural = pd.read_parquet(
        ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet"
    )
    natural_sequences = (
        natural[["protein_id", "protein_sequence"]]
        .drop_duplicates("protein_id")
        .set_index("protein_id")["protein_sequence"]
        .to_dict()
    )
    designed_sequences = (
        df[["protein_id", "protein_sequence"]]
        .drop_duplicates("protein_id")
        .set_index("protein_id")["protein_sequence"]
        .to_dict()
    )
    splits = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "natural_pbm_splits.csv")
    train_ids = splits.loc[splits["split"].eq("train"), "protein_id"].tolist()
    train_ids = [protein for protein in train_ids if protein in natural_sequences]

    embeddings_path = (
        ROOT
        / "data"
        / "interim"
        / "v0_4_2"
        / "frozen_plm_embeddings_esm2_t12_35M_UR50D.parquet"
    )
    embeddings = pd.read_parquet(embeddings_path)
    emb_cols = [column for column in embeddings.columns if column.startswith("emb_")]
    emb_lookup = {
        row["protein_id"]: row[emb_cols].to_numpy(dtype=float)
        for _, row in embeddings.iterrows()
    }

    rows = []
    for protein_id, group in df.groupby("protein_id", sort=True):
        designed_seq = designed_sequences[protein_id]
        identities = [
            sequence_identity(designed_seq, natural_sequences[train_id])
            for train_id in train_ids
        ]
        nearest_embedding = [
            float(np.linalg.norm(emb_lookup[protein_id] - emb_lookup[train_id]))
            for train_id in train_ids
            if protein_id in emb_lookup and train_id in emb_lookup
        ]
        top5 = group.nlargest(max(1, int(math.ceil(len(group) * 0.05)),), "experimental_score")
        motif = str(group["designed_binding_site_motif"].iloc[0])
        motif_len = len(motif)
        motif_distances = group["motif_hamming_distance"].to_numpy(dtype=float)
        rows.append(
            {
                "protein_id": protein_id,
                "motif_sequence": motif,
                "motif_length": motif_len,
                "protein_sequence_cluster": groups.loc[
                    groups["protein_id"].eq(protein_id), "protein_sequence_cluster"
                ].iloc[0],
                "original_target_group": groups.loc[
                    groups["protein_id"].eq(protein_id), "original_target_group"
                ].iloc[0],
                "assay_target_group": groups.loc[
                    groups["protein_id"].eq(protein_id), "assay_target_group"
                ].iloc[0],
                "motif_group": groups.loc[
                    groups["protein_id"].eq(protein_id), "motif_group"
                ].iloc[0],
                "max_natural_train_sequence_identity": max(identities) if identities else np.nan,
                "nearest_natural_train_esm_euclidean_distance": (
                    min(nearest_embedding) if nearest_embedding else np.nan
                ),
                "experimental_score_std": score_distribution_features(group)["experimental_score_std"],
                "experimental_score_variance": score_distribution_features(group)["experimental_score_variance"],
                "experimental_score_iqr": score_distribution_features(group)["experimental_score_iqr"],
                "experimental_score_q95": score_distribution_features(group)["experimental_score_q95"],
                "experimental_score_q99": score_distribution_features(group)["experimental_score_q99"],
                "top_1pct_mean": score_distribution_features(group)["top_1pct_mean"],
                "top_1pct_minus_median": score_distribution_features(group)["top_1pct_minus_median"],
                "top_1pct_shifted_mass_fraction": score_distribution_features(group)["top_1pct_shifted_mass_fraction"],
                "score_histogram_entropy_normalized": score_distribution_features(group)["score_histogram_entropy_normalized"],
                "n_unique_experimental_scores": score_distribution_features(group)["n_unique_experimental_scores"],
                "top5_mean_motif_hamming_distance": float(np.mean(top5["motif_hamming_distance"])),
                "top5_fraction_within_one_motif_mismatch": float(
                    (top5["motif_hamming_distance"] <= 1).mean()
                ),
                "all_rows_mean_motif_hamming_distance": float(np.mean(motif_distances)),
                "all_rows_median_motif_hamming_distance": float(np.median(motif_distances)),
            }
        )
    out = pd.DataFrame(rows)
    out = out.merge(seq_perf, on="protein_id", how="left")
    out = out.merge(pc_perf, on="protein_id", how="left")
    out = out.merge(plm_perf, on="protein_id", how="left")
    out.to_csv(V042_TABLES / "designed_difficulty_factors.csv", index=False)
    return out


def write_figures(df: pd.DataFrame, common_hard: pd.DataFrame, distance: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    factors = pd.read_csv(V042_TABLES / "designed_difficulty_factors.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.scatter(
        factors["max_natural_train_sequence_identity"],
        factors["frozen_plm_spearman"],
        color="#4C78A8",
        s=55,
    )
    for _, row in factors.iterrows():
        ax.annotate(row["protein_id"], (row["max_natural_train_sequence_identity"], row["frozen_plm_spearman"]), fontsize=8, xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("Maximum identity to natural train protein")
    ax.set_ylabel("FrozenPLM designed Spearman")
    ax.set_title("Designed performance vs natural-train sequence similarity")
    fig.tight_layout()
    fig.savefig(V042_FIGURES / "fig_v0_4_2_6_performance_vs_train_similarity.png", dpi=300)
    plt.close(fig)

    plotted = distance[distance["method"].isin(["sequence_kmer3", "SimpleProteinConditionalBaseline", "FrozenPLMProteinConditionalBaseline"])]
    summary = plotted.groupby(["method", "motif_distance_bin"], as_index=False)["spearman"].median()
    order = ["0", "1", "2", "3+"]
    x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    colors = {
        "sequence_kmer3": "#4C78A8",
        "SimpleProteinConditionalBaseline": "#54A24B",
        "FrozenPLMProteinConditionalBaseline": "#F58518",
    }
    for method, group in summary.groupby("method"):
        values = [
            group.loc[group["motif_distance_bin"].eq(level), "spearman"].iloc[0]
            if (group["motif_distance_bin"] == level).any()
            else np.nan
            for level in order
        ]
        ax.plot(x, values, marker="o", label=method, color=colors[method])
    ax.set_xticks(x, order)
    ax.set_xlabel("RC-aware motif Hamming distance bin")
    ax.set_ylabel("Median per-protein Spearman")
    ax.set_title("Ranking performance by motif-distance regime")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(V042_FIGURES / "fig_v0_4_2_7_performance_by_motif_distance.png", dpi=300)
    plt.close(fig)

    counts = (
        common_hard.query("failure_type == 'high_experiment_low_all_core'")
        .groupby("protein_id")
        .size()
        .reindex(DESIGNED_PROTEINS, fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.bar(counts.index, counts.values, color="#D45087")
    ax.set_xlabel("Designed DBP")
    ax.set_ylabel("Number of common hard cases")
    ax.set_title("Core common-hard experimental specificity cases")
    fig.tight_layout()
    fig.savefig(V042_FIGURES / "fig_v0_4_2_8_common_hard_set.png", dpi=300)
    plt.close(fig)


def main() -> None:
    designed = load_designed()
    scored = add_descriptive_features(load_predictions(designed))
    resolution = write_disagreement_resolution(scored)
    common_hard = write_common_hard_set(scored)
    distance = write_performance_by_distance(scored)
    factors = write_difficulty_factors(scored)

    summary = {
        "n_designed_rows": int(len(scored)),
        "n_disagreement_candidates": int(disagreement_mask(scored).sum()),
        "n_common_high_experiment_low_all_core": int(
            (common_hard["failure_type"] == "high_experiment_low_all_core").sum()
        ),
        "n_common_low_experiment_high_all_core": int(
            (common_hard["failure_type"] == "low_experiment_high_all_core").sum()
        ),
        "n_protein_sequence_clusters": int(
            factors["protein_sequence_cluster"].nunique()
        ),
        "n_original_target_groups": int(factors["original_target_group"].nunique()),
        "n_assay_target_groups": int(factors["assay_target_group"].nunique()),
        "n_motif_groups": int(factors["motif_group"].nunique()),
    }
    pd.DataFrame([summary]).to_json(
        V042_TABLES / "v0_4_2_diagnostic_summary.json", orient="records", indent=2
    )
    write_figures(scored, common_hard, distance)
    print(pd.DataFrame([summary]).to_string(index=False))
    print("\nDisagreement resolution:\n", resolution.to_string(index=False))
    print("\nDifficulty factors:\n", factors.to_string(index=False))


if __name__ == "__main__":
    main()
