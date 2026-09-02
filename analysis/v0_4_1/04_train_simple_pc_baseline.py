from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.models.simple_protein_conditional_baseline import dna_kmer_features, protein_composition_features
from src.utils import ensure_dir, project_root
from src.v0_4_evaluation import compute_ranking_metrics


ROOT = project_root()
PROCESSED = ROOT / "data" / "processed"
RESULTS = ensure_dir(ROOT / "results" / "v0_4_1")
TABLES = ensure_dir(RESULTS / "tables")
FIGURES = ensure_dir(RESULTS / "figures")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4_1")
DOCS = ensure_dir(ROOT / "docs" / "v0_4_1")
RNG = np.random.default_rng(42)
ALPHAS = [0.1, 1.0, 10.0, 100.0]
TRAIN_ROWS_PER_PROTEIN = 2500
VALID_ROWS_PER_PROTEIN = 4000


def feature_vector(protein_sequence: str, dna_sequence: str) -> np.ndarray:
    p = protein_composition_features(protein_sequence)
    d = dna_kmer_features(dna_sequence)
    p_small = p[[20, 21, 22, 23, 24]]
    return np.concatenate([p, d, np.outer(p_small, d).ravel()]).astype(np.float32)


def build_feature_matrix(df: pd.DataFrame, seqs: dict[str, str], dna_col: str) -> np.ndarray:
    rows = []
    for protein_id, dna in zip(df["protein_id"], df[dna_col]):
        rows.append(feature_vector(seqs[protein_id], str(dna)))
    return np.vstack(rows)


def stratified_sample(df: pd.DataFrame, n_per_protein: int) -> pd.DataFrame:
    parts = []
    for _, group in df.groupby("protein_id", sort=True):
        n = min(n_per_protein, len(group))
        parts.append(group.sample(n=n, random_state=42))
    return pd.concat(parts, ignore_index=True)


def fit_ridge(x: np.ndarray, y: np.ndarray, alpha: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    mu = x.mean(axis=0)
    sigma = x.std(axis=0)
    sigma[sigma == 0] = 1.0
    xs = (x - mu) / sigma
    y_mean = float(y.mean())
    yc = y - y_mean
    xtx = xs.T @ xs
    xtx.flat[:: xtx.shape[0] + 1] += alpha
    w = np.linalg.solve(xtx, xs.T @ yc)
    return w.astype(np.float32), mu.astype(np.float32), sigma.astype(np.float32), y_mean


def predict(df: pd.DataFrame, seqs: dict[str, str], dna_col: str, model: dict, chunk: int = 50_000) -> np.ndarray:
    preds = []
    for start in range(0, len(df), chunk):
        part = df.iloc[start : start + chunk]
        x = build_feature_matrix(part, seqs, dna_col)
        xs = (x - model["mu"]) / model["sigma"]
        preds.append(xs @ model["weights"] + model["y_mean"])
    return np.concatenate(preds)


def evaluate(df: pd.DataFrame, score_col: str, pred_col: str, dataset: str) -> pd.DataFrame:
    rows = []
    for protein_id, group in df.groupby("protein_id", sort=True):
        eval_df = group.rename(columns={"canonical_rc": "canonical_7mer"})
        metrics = compute_ranking_metrics(eval_df, score_col, pred_col)
        rows.append(
            {
                "dataset": dataset,
                "protein_id": protein_id,
                "baseline": "SimpleProteinConditionalBaseline_composition_ridge",
                "spearman": metrics.spearman,
                "ndcg_1pct": metrics.ndcg_1pct,
                "ndcg_5pct": metrics.ndcg_5pct,
                "pairwise_accuracy": metrics.pairwise_accuracy,
                "top1pct_recovery": metrics.top1pct_recovery,
                "n_rc_classes": metrics.n_rc_classes,
            }
        )
    return pd.DataFrame(rows)


def macro_by_dataset(perf: pd.DataFrame) -> pd.DataFrame:
    metric_cols = ["spearman", "ndcg_1pct", "ndcg_5pct", "pairwise_accuracy", "top1pct_recovery"]
    rows = []
    for (dataset, baseline), group in perf.groupby(["dataset", "baseline"], sort=True):
        for metric in metric_cols:
            vals = group[metric].dropna().to_numpy(dtype=float)
            rows.append(
                {
                    "dataset": dataset,
                    "baseline": baseline,
                    "metric": metric,
                    "n_proteins": int(group["protein_id"].nunique()),
                    "n_proteins_with_metric": int(vals.size),
                    "mean": float(np.mean(vals)) if vals.size else np.nan,
                    "median": float(np.median(vals)) if vals.size else np.nan,
                    "std": float(np.std(vals, ddof=1)) if vals.size > 1 else np.nan,
                }
            )
    return pd.DataFrame(rows)


def load_designed() -> pd.DataFrame:
    designed = pd.read_parquet(PROCESSED / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet")
    designed = designed.rename(
        columns={
            "canonical_7mer": "canonical_rc",
            "experimental_escore_consensus": "experimental_score",
        }
    )
    designed["experimental_percentile"] = designed.groupby("protein_id")["experimental_score"].rank(pct=True, ascending=True)
    return designed


def disagreement_resolution(designed_pred: pd.DataFrame) -> pd.DataFrame:
    scored = pd.read_parquet(PROCESSED / "v0_3_1" / "designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1.parquet")
    counts = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "all_disagreement_candidate_counts.csv")
    pred = designed_pred[["protein_id", "canonical_rc", "simple_pc_score"]].rename(columns={"canonical_rc": "canonical_7mer"})
    merged = scored.merge(pred, on=["protein_id", "canonical_7mer"], how="left")
    rows = []
    examples = []
    for _, c in counts.iterrows():
        protein_id = c["protein_id"]
        group = merged[merged["protein_id"] == protein_id].copy()
        group["simple_pc_percentile"] = group["simple_pc_score"].rank(pct=True, ascending=True)
        is_disagreement = (
            (group["experimental_escore_consensus"] >= float(c["experimental_score_threshold"]) - 1e-12)
            & (group["hamming_similarity_to_paper_motif_rc_aware"] <= float(c["sequence_similarity_threshold"]) + 1e-12)
        )
        cand = group[is_disagreement].copy()
        resolved = cand["simple_pc_percentile"] >= 0.90
        rows.append(
            {
                "protein_id": protein_id,
                "baseline": "SimpleProteinConditionalBaseline_composition_ridge",
                "n_v0_3_1_disagreement_candidates": int(c["n_disagreement"]),
                "n_evaluable_candidates": int(cand["simple_pc_score"].notna().sum()),
                "n_resolved": int(resolved.sum()),
                "resolution_rate": float(resolved.mean()) if len(cand) else np.nan,
                "resolution_definition": "candidate is sequence-vs-experiment disagreement and SimplePC ranks it in top 10% within protein",
            }
        )
        top = cand.assign(simple_pc_resolved=resolved).sort_values("experimental_escore_consensus", ascending=False).head(20)
        examples.append(top)
    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "simple_pc_disagreement_resolution.csv", index=False)
    if examples:
        pd.concat(examples, ignore_index=True).to_csv(TABLES / "simple_pc_disagreement_examples.csv", index=False)
    return out


def make_figures(perf: pd.DataFrame, natural_pred: pd.DataFrame, designed_pred: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    macro = macro_by_dataset(perf)
    spearman = macro[macro["metric"] == "spearman"].copy()
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.bar(spearman["baseline"] + "\n" + spearman["dataset"], spearman["median"], color=["#3A6EA5", "#C95D63"])
    ax.axhline(0.591, color="black", linestyle="--", linewidth=1, label="Designed uPBM replicate Spearman reference")
    ax.set_ylabel("Macro median Spearman")
    ax.set_title("Simple PC Natural Held-out vs Designed External")
    ax.tick_params(axis="x", labelrotation=20)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_v0_4_1_4_simple_pc_natural_vs_designed.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    natural_example = natural_pred[natural_pred["protein_id"] == sorted(natural_pred["protein_id"].unique())[0]].sample(
        n=min(4000, len(natural_pred[natural_pred["protein_id"] == sorted(natural_pred["protein_id"].unique())[0]])),
        random_state=42,
    )
    designed_example = designed_pred[designed_pred["protein_id"] == "DBP1"].sample(n=4000, random_state=42)
    axes[0].scatter(natural_example["simple_pc_score"], natural_example["experimental_score"], s=4, alpha=0.35)
    axes[0].set_title("Natural held-out example")
    axes[0].set_xlabel("Simple PC score")
    axes[0].set_ylabel("PBM E-score")
    axes[1].scatter(designed_example["simple_pc_score"], designed_example["experimental_score"], s=4, alpha=0.35, color="#C95D63")
    axes[1].set_title("Designed external example")
    axes[1].set_xlabel("Simple PC score")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_v0_4_1_7_failure_case_comparison.png", dpi=300)
    plt.close(fig)


def main() -> None:
    natural = pd.read_parquet(PROCESSED / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    splits = pd.read_csv(METADATA / "natural_pbm_splits.csv")
    seqs = dict(natural[["protein_id", "protein_sequence"]].drop_duplicates().itertuples(index=False, name=None))
    natural = natural.merge(splits[["protein_id", "split"]], on="protein_id", how="left")
    train = stratified_sample(natural[natural["split"] == "train"], TRAIN_ROWS_PER_PROTEIN)
    valid = stratified_sample(natural[natural["split"] == "validation"], VALID_ROWS_PER_PROTEIN)
    x_train = build_feature_matrix(train, seqs, "canonical_rc")
    y_train = train["experimental_percentile"].to_numpy(dtype=np.float32)
    x_valid = build_feature_matrix(valid, seqs, "canonical_rc")
    y_valid = valid["experimental_percentile"].to_numpy(dtype=np.float32)

    records = []
    best_model = None
    best_alpha = None
    best_spearman = -np.inf
    for alpha in ALPHAS:
        weights, mu, sigma, y_mean = fit_ridge(x_train, y_train, alpha)
        model = {"weights": weights, "mu": mu, "sigma": sigma, "y_mean": y_mean}
        pred = ((x_valid - mu) / sigma) @ weights + y_mean
        val_eval = valid.assign(simple_pc_score=pred)
        perf = evaluate(val_eval, "experimental_score", "simple_pc_score", "validation")
        median_spearman = float(perf["spearman"].median())
        records.append({"alpha": alpha, "validation_macro_median_spearman": median_spearman})
        if median_spearman > best_spearman:
            best_spearman = median_spearman
            best_alpha = alpha
            best_model = model
    assert best_model is not None

    natural_test = natural[natural["split"] == "natural_test"].copy()
    natural_test["simple_pc_score"] = predict(natural_test, seqs, "canonical_rc", best_model)
    designed = load_designed()
    designed_seqs = dict(designed[["protein_id", "protein_sequence"]].drop_duplicates().itertuples(index=False, name=None))
    designed["simple_pc_score"] = predict(designed, designed_seqs, "canonical_rc", best_model)

    natural_test.to_parquet(TABLES / "simple_pc_natural_test_predictions.parquet", index=False)
    designed.to_parquet(TABLES / "simple_pc_designed_predictions.parquet", index=False)
    natural_perf = evaluate(natural_test, "experimental_score", "simple_pc_score", "natural_test")
    designed_perf = evaluate(designed, "experimental_score", "simple_pc_score", "designed_external")
    perf = pd.concat([natural_perf, designed_perf], ignore_index=True)
    perf.to_csv(TABLES / "simple_pc_performance.csv", index=False)
    macro = macro_by_dataset(perf)
    macro.to_csv(TABLES / "simple_pc_performance_macro.csv", index=False)
    pd.DataFrame(records).to_csv(TABLES / "simple_pc_alpha_selection.csv", index=False)
    disagreement_resolution(designed)
    make_figures(perf, natural_test, designed)
    model_meta = {
        "baseline": "SimpleProteinConditionalBaseline_composition_ridge",
        "status": "trained",
        "not_proposed_method": True,
        "protein_representation": "amino-acid composition and simple physicochemical summaries; no protein LM fine-tuning",
        "dna_representation": "padded one-hot, base composition, dinucleotide composition, GC, DNA length",
        "interaction": "outer product of five protein summary features with DNA features",
        "training_data": "UniPROBE natural PBM train proteins only",
        "target": "within-protein percentile of processed 8-mer PBM E-score",
        "selected_alpha": best_alpha,
        "validation_macro_median_spearman": best_spearman,
        "train_rows_per_protein": TRAIN_ROWS_PER_PROTEIN,
        "validation_rows_per_protein": VALID_ROWS_PER_PROTEIN,
        "seed": 42,
    }
    (RESULTS / "simple_pc_model_metadata.json").write_text(json.dumps(model_meta, indent=2), encoding="utf-8")
    print(json.dumps(model_meta, indent=2))
    print(macro.to_string(index=False))


if __name__ == "__main__":
    main()
