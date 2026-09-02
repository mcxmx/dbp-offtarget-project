from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.models.simple_protein_conditional_baseline import dna_kmer_features
from src.utils import ensure_dir, project_root
from src.v0_4_evaluation import compute_ranking_metrics


ROOT = project_root()
TODAY = date.today().isoformat()
PROCESSED = ROOT / "data" / "processed"
V042_DATA = ensure_dir(PROCESSED / "v0_4_2")
V042_INTERIM = ensure_dir(ROOT / "data" / "interim" / "v0_4_2")
V042_META = ensure_dir(ROOT / "metadata" / "v0_4_2")
V042_RESULTS = ensure_dir(ROOT / "results" / "v0_4_2")
V042_TABLES = ensure_dir(V042_RESULTS / "tables")
V042_FIGURES = ensure_dir(V042_RESULTS / "figures")
V042_DOCS = ensure_dir(ROOT / "docs" / "v0_4_2")

SEED = 42
CHECKPOINT_NAME = "esm2_t12_35M_UR50D"
ESM_REPR_LAYER = 12
MAX_ESM_RESIDUES = 1022
TRAIN_ROWS_PER_PROTEIN = 1800
VALID_ROWS_PER_PROTEIN = 4000
PREDICT_CHUNK = 50_000
ALPHAS = [1.0, 10.0, 100.0, 1000.0]
AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def clean_for_esm(seq: str) -> tuple[str, int, int]:
    raw = str(seq).upper()
    replaced = sum(aa not in AA20 for aa in raw)
    cleaned = "".join(aa if aa in AA20 else "X" for aa in raw)
    truncated = max(0, len(cleaned) - MAX_ESM_RESIDUES)
    return cleaned[:MAX_ESM_RESIDUES], replaced, truncated


def load_designed() -> pd.DataFrame:
    designed = pd.read_parquet(PROCESSED / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet")
    designed = designed.rename(
        columns={
            "canonical_7mer": "canonical_rc",
            "experimental_escore_consensus": "experimental_score",
        }
    )
    designed["experimental_percentile"] = designed.groupby("protein_id")["experimental_score"].rank(
        pct=True, ascending=True
    )
    designed["dna_length"] = 7
    designed["dataset"] = "designed_external"
    return designed


def unique_protein_sequences(natural: pd.DataFrame, designed: pd.DataFrame) -> pd.DataFrame:
    nat = natural[["protein_id", "protein_sequence"]].drop_duplicates().assign(sequence_dataset="natural_pbm")
    des = designed[["protein_id", "protein_sequence"]].drop_duplicates().assign(sequence_dataset="designed_upbm")
    seqs = pd.concat([nat, des], ignore_index=True)
    return seqs.drop_duplicates("protein_id")


def compute_esm_embeddings(seqs: pd.DataFrame) -> pd.DataFrame:
    out_path = V042_INTERIM / f"frozen_plm_embeddings_{CHECKPOINT_NAME}.parquet"
    manifest_path = V042_META / f"frozen_plm_embedding_manifest_{CHECKPOINT_NAME}.csv"
    if out_path.exists() and manifest_path.exists():
        return pd.read_parquet(out_path)

    import esm

    torch.manual_seed(SEED)
    model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    rows = []
    manifest_rows = []
    with torch.no_grad():
        for _, row in seqs.sort_values("protein_id").iterrows():
            cleaned, replaced, truncated = clean_for_esm(row["protein_sequence"])
            labels, strings, tokens = batch_converter([(row["protein_id"], cleaned)])
            result = model(tokens, repr_layers=[ESM_REPR_LAYER], return_contacts=False)
            reps = result["representations"][ESM_REPR_LAYER][0, 1 : len(cleaned) + 1]
            emb = reps.mean(dim=0).cpu().numpy().astype(np.float32)
            rows.append(
                {
                    "protein_id": row["protein_id"],
                    "sequence_dataset": row["sequence_dataset"],
                    **{f"emb_{i:03d}": float(value) for i, value in enumerate(emb)},
                }
            )
            manifest_rows.append(
                {
                    "protein_id": row["protein_id"],
                    "checkpoint_name": CHECKPOINT_NAME,
                    "repr_layer": ESM_REPR_LAYER,
                    "embedding_dim": int(emb.shape[0]),
                    "raw_sequence_length": int(len(str(row["protein_sequence"]))),
                    "embedded_sequence_length": int(len(cleaned)),
                    "n_noncanonical_replaced_with_X": int(replaced),
                    "n_truncated_residues": int(truncated),
                    "protein_lm_frozen": True,
                    "notes": "Mean-pooled residue representation; protein LM not fine-tuned.",
                }
            )
    emb_df = pd.DataFrame(rows)
    emb_df.to_parquet(out_path, index=False)
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)
    return emb_df


def sampled_rows(df: pd.DataFrame, n_per_protein: int) -> pd.DataFrame:
    parts = []
    for _, group in df.groupby("protein_id", sort=True):
        n = min(n_per_protein, len(group))
        parts.append(group.sample(n=n, random_state=SEED))
    return pd.concat(parts, ignore_index=True)


def embedding_lookup(embeddings: pd.DataFrame) -> dict[str, np.ndarray]:
    emb_cols = [c for c in embeddings.columns if c.startswith("emb_")]
    return {
        str(row["protein_id"]): row[emb_cols].to_numpy(dtype=np.float32)
        for _, row in embeddings.iterrows()
    }


def feature_matrix(df: pd.DataFrame, emb: dict[str, np.ndarray], dna_col: str = "canonical_rc") -> np.ndarray:
    rows = []
    for protein_id, dna in zip(df["protein_id"], df[dna_col]):
        p = emb[str(protein_id)]
        d = dna_kmer_features(str(dna))
        rows.append(np.concatenate([p, d]))
    return np.vstack(rows).astype(np.float32)


def fit_ridge(x: np.ndarray, y: np.ndarray, alpha: float) -> dict[str, np.ndarray | float]:
    mu = x.mean(axis=0)
    sigma = x.std(axis=0)
    sigma[sigma == 0] = 1.0
    xs = (x - mu) / sigma
    y_mean = float(y.mean())
    yc = y - y_mean
    xtx = xs.T @ xs
    xtx.flat[:: xtx.shape[0] + 1] += alpha
    weights = np.linalg.solve(xtx, xs.T @ yc)
    return {
        "weights": weights.astype(np.float32),
        "feature_mean": mu.astype(np.float32),
        "feature_scale": sigma.astype(np.float32),
        "intercept": y_mean,
        "alpha": alpha,
    }


def predict(df: pd.DataFrame, emb: dict[str, np.ndarray], model: dict) -> np.ndarray:
    preds = []
    for start in range(0, len(df), PREDICT_CHUNK):
        chunk = df.iloc[start : start + PREDICT_CHUNK]
        x = feature_matrix(chunk, emb)
        xs = (x - model["feature_mean"]) / model["feature_scale"]
        preds.append(xs @ model["weights"] + model["intercept"])
    return np.concatenate(preds).astype(np.float32)


def evaluate(df: pd.DataFrame, pred_col: str, dataset: str) -> pd.DataFrame:
    rows = []
    for protein_id, group in df.groupby("protein_id", sort=True):
        eval_df = group.rename(columns={"canonical_rc": "canonical_7mer"})
        metrics = compute_ranking_metrics(eval_df, "experimental_score", pred_col)
        rows.append(
            {
                "dataset": dataset,
                "protein_id": protein_id,
                "baseline": "FrozenPLMProteinConditionalBaseline_ESM2_t12_ridge",
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
    rows = []
    for (dataset, baseline), group in perf.groupby(["dataset", "baseline"], sort=True):
        for metric in ["spearman", "ndcg_1pct", "ndcg_5pct", "pairwise_accuracy", "top1pct_recovery"]:
            values = group[metric].dropna().to_numpy(dtype=float)
            rows.append(
                {
                    "dataset": dataset,
                    "baseline": baseline,
                    "metric": metric,
                    "n_proteins": int(group["protein_id"].nunique()),
                    "n_proteins_with_metric": int(values.size),
                    "mean": float(np.mean(values)) if values.size else np.nan,
                    "median": float(np.median(values)) if values.size else np.nan,
                    "std": float(np.std(values, ddof=1)) if values.size > 1 else np.nan,
                }
            )
    return pd.DataFrame(rows)


def write_model_artifacts(model: dict, validation_record: dict) -> None:
    np.savez_compressed(
        V042_DATA / "frozen_plm_primary_head_v0_4_2.npz",
        weights=model["weights"],
        feature_mean=model["feature_mean"],
        feature_scale=model["feature_scale"],
        intercept=np.array([model["intercept"]], dtype=np.float32),
    )
    metadata = {
        "baseline": "FrozenPLMProteinConditionalBaseline",
        "not_proposed_method": True,
        "protein_encoder": CHECKPOINT_NAME,
        "protein_encoder_frozen": True,
        "repr_layer": ESM_REPR_LAYER,
        "protein_representation": "mean-pooled frozen ESM-2 residue embeddings",
        "protein_sequence_version": "FULL_LENGTH_REFERENCE for natural PBM; designed DBP sequences from v0.3.1",
        "construct_aligned_training": False,
        "dna_representation": "one-hot and simple k-mer/composition features for 7-mer/8-mer DNA",
        "interaction_head": "raw mean-pooled ESM embedding concatenated with DNA k-mer features; ridge regression",
        "training_data": "natural PBM train split only",
        "validation_data": "natural PBM validation split only",
        "external_test": "GSE237017 designed uPBM, not used for alpha selection",
        "target": "within-protein experimental percentile of processed PBM E-score",
        "selected_alpha": validation_record["alpha"],
        "validation_macro_median_spearman": validation_record["validation_macro_median_spearman"],
        "seed": SEED,
        "date": TODAY,
    }
    (V042_RESULTS / "frozen_plm_model_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def make_figures(perf: pd.DataFrame, natural_pred: pd.DataFrame, designed_pred: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    macro = macro_by_dataset(perf)
    spearman = macro[macro["metric"].eq("spearman")].copy()
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    labels = spearman["dataset"].str.replace("_", " ", regex=False)
    ax.bar(labels, spearman["median"], color=["#4C78A8", "#F58518"])
    ax.axhline(0.591, color="black", linestyle="--", linewidth=1, label="Designed uPBM replicate reference")
    ax.set_ylabel("Macro median Spearman")
    ax.set_title("Frozen ESM-2 Baseline")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(V042_FIGURES / "fig_v0_4_2_2_natural_vs_designed_frozen_plm.png", dpi=300)
    plt.close(fig)

    designed_perf = perf[perf["dataset"].eq("designed_external")].copy()
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.bar(designed_perf["protein_id"], designed_perf["spearman"], color="#F58518")
    ax.axhline(0.232, color="#4C78A8", linestyle=":", linewidth=1.2, label="Best sequence-only median")
    ax.axhline(0.362, color="#54A24B", linestyle="-.", linewidth=1.2, label="SimplePC median")
    ax.axhline(0.591, color="black", linestyle="--", linewidth=1.0, label="Replicate reference")
    ax.set_ylabel("Spearman")
    ax.set_title("Designed DBP FrozenPLM Performance")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(V042_FIGURES / "fig_v0_4_2_4_per_designed_dbp_frozen_plm.png", dpi=300)
    plt.close(fig)

    for protein_id in ["DBP6", "DBP48"]:
        group = designed_pred[designed_pred["protein_id"].eq(protein_id)].sample(
            n=min(4000, int(designed_pred["protein_id"].eq(protein_id).sum())),
            random_state=SEED,
        )
        fig, ax = plt.subplots(figsize=(4.5, 4.0))
        ax.scatter(group["frozen_plm_score"], group["experimental_score"], s=5, alpha=0.35, color="#E45756")
        ax.set_xlabel("FrozenPLM score")
        ax.set_ylabel("processed uPBM E-score")
        ax.set_title(f"{protein_id}: experimental vs FrozenPLM")
        fig.tight_layout()
        fig.savefig(V042_FIGURES / f"fig_v0_4_2_5_{protein_id}_experimental_vs_frozen_plm.png", dpi=300)
        plt.close(fig)


def main() -> None:
    natural = pd.read_parquet(PROCESSED / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    splits = pd.read_csv(ROOT / "metadata" / "v0_4_1" / "natural_pbm_splits.csv")
    natural = natural.merge(splits[["protein_id", "split"]], on="protein_id", how="left")
    designed = load_designed()
    proteins = unique_protein_sequences(natural, designed)
    embeddings = compute_esm_embeddings(proteins)
    emb = embedding_lookup(embeddings)

    train = sampled_rows(natural[natural["split"].eq("train")], TRAIN_ROWS_PER_PROTEIN)
    valid = sampled_rows(natural[natural["split"].eq("validation")], VALID_ROWS_PER_PROTEIN)
    x_train = feature_matrix(train, emb)
    y_train = train["experimental_percentile"].to_numpy(dtype=np.float32)
    x_valid = feature_matrix(valid, emb)

    records = []
    best_model: dict | None = None
    best_record: dict | None = None
    for alpha in ALPHAS:
        model = fit_ridge(x_train, y_train, alpha)
        valid_pred = ((x_valid - model["feature_mean"]) / model["feature_scale"]) @ model["weights"] + model["intercept"]
        valid_eval = valid.assign(frozen_plm_score=valid_pred)
        valid_perf = evaluate(valid_eval, "frozen_plm_score", "validation")
        median_spearman = float(valid_perf["spearman"].median())
        record = {"alpha": alpha, "validation_macro_median_spearman": median_spearman}
        records.append(record)
        if best_record is None or median_spearman > best_record["validation_macro_median_spearman"]:
            best_record = record
            best_model = model
    assert best_model is not None and best_record is not None

    pd.DataFrame(records).to_csv(V042_TABLES / "frozen_plm_alpha_selection.csv", index=False)
    write_model_artifacts(best_model, best_record)

    natural_test = natural[natural["split"].eq("natural_test")].copy()
    natural_test["frozen_plm_score"] = predict(natural_test, emb, best_model)
    designed["frozen_plm_score"] = predict(designed, emb, best_model)
    natural_test.to_parquet(V042_TABLES / "frozen_plm_natural_test_predictions.parquet", index=False)
    designed.to_parquet(V042_TABLES / "frozen_plm_designed_predictions.parquet", index=False)

    perf = pd.concat(
        [
            evaluate(natural_test, "frozen_plm_score", "natural_test"),
            evaluate(designed, "frozen_plm_score", "designed_external"),
        ],
        ignore_index=True,
    )
    perf.to_csv(V042_TABLES / "frozen_plm_performance.csv", index=False)
    macro = macro_by_dataset(perf)
    macro.to_csv(V042_TABLES / "frozen_plm_performance_macro.csv", index=False)
    make_figures(perf, natural_test, designed)

    (V042_DOCS / "FROZEN_PLM_BASELINE.md").write_text(
        f"""# v0.4.2 FrozenPLMProteinConditionalBaseline

Date: {TODAY}

This baseline uses frozen ESM-2 `{CHECKPOINT_NAME}` mean-pooled protein embeddings concatenated with simple DNA features and scored by ridge regression. It is a baseline only, not the proposed model.

## Training and validation

- Training data: natural PBM train proteins only.
- Validation: natural PBM validation proteins only.
- External test: GSE237017 designed DBPs only after alpha selection.
- Selected alpha: {best_record['alpha']}
- Validation macro median Spearman: {best_record['validation_macro_median_spearman']:.3f}

The protein LM is not fine-tuned. Designed DBP rows do not influence checkpoint, alpha, feature encoding, or any hyperparameter selection.
""",
        encoding="utf-8",
    )
    print(json.dumps({"selected_alpha": best_record["alpha"], "validation": best_record}, indent=2))
    print(macro.to_string(index=False))


if __name__ == "__main__":
    main()
