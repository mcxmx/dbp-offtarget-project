from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from torch import nn

from src.v0_4_evaluation import compute_ranking_metrics
from src.v0_5_models import (
    CandidateDNAOnly,
    ProteinCandidate,
    ProteinTargetCandidate,
    TargetCandidateOnly,
    batch_rc_symmetric_one_hot,
    model_parameter_counts,
    score_rc,
    target_edit_control,
    target_hamming_control,
    target_kmer_overlap_control,
    target_window_features,
)
from src.utils import ensure_dir, project_root


ROOT = project_root()
DESIGNED_PATH = ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet"
TARGET_PATH = ROOT / "metadata" / "v0_5" / "designed_target_manifest_v0_5.csv"
EMBEDDING_PATH = ROOT / "data" / "interim" / "v0_4_2" / "frozen_plm_embeddings_esm2_t12_35M_UR50D.parquet"
SPLIT_PATH = ROOT / "metadata" / "v0_5" / "v0_5_split_manifest.csv"


@dataclass(frozen=True)
class V05Config:
    seed: int = 42
    evaluation_seeds: tuple[int, ...] = (17, 29, 43)
    smoke_split_name: str = "protein_cluster_loco"
    smoke_fold_id: str = "protein_cluster_loco_fold_1"
    hidden_dim: int = 32
    capacity_matched_hidden_dim: int = 40
    learning_rate: float = 0.01
    weight_decay: float = 1e-4
    epochs: int = 18
    pair_count_per_protein: int = 512
    batch_size: int = 128
    tie_tolerance: float = 1e-10

    def as_dict(self) -> dict[str, Any]:
        return {
            "seed": self.seed,
            "evaluation_seeds": list(self.evaluation_seeds),
            "smoke_split_name": self.smoke_split_name,
            "smoke_fold_id": self.smoke_fold_id,
            "hidden_dim": self.hidden_dim,
            "capacity_matched_hidden_dim": self.capacity_matched_hidden_dim,
            "learning_rate": self.learning_rate,
            "weight_decay": self.weight_decay,
            "epochs": self.epochs,
            "pair_count_per_protein": self.pair_count_per_protein,
            "batch_size": self.batch_size,
            "tie_tolerance": self.tie_tolerance,
            "protein_representation": "frozen ESM-2 esm2_t12_35M_UR50D mean-pooled embeddings",
            "protein_embedding_dim": 480,
            "dna_representation": "RC-symmetric one-hot 7-mer features",
            "target_representation": "all canonical RC target windows of length 7, mean pooled",
            "primary_loss": "within-protein logistic pairwise ranking loss",
            "pair_sampling": "deterministic 40% easy, 35% medium, 25% hard by within-protein rank difference",
            "designed_test_used_for_selection": False,
            "status": "smoke_training_only_not_primary_scientific_result",
        }


def set_seed(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)


def load_v05_data() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray], pd.DataFrame]:
    benchmark = pd.read_parquet(DESIGNED_PATH).rename(
        columns={
            "canonical_7mer": "candidate_dna",
            "experimental_escore_consensus": "experimental_score",
        }
    )
    benchmark = benchmark[["protein_id", "candidate_dna", "experimental_score"]].copy()
    targets = pd.read_csv(TARGET_PATH)
    embeddings = pd.read_parquet(EMBEDDING_PATH)
    embedding_columns = [column for column in embeddings.columns if column.startswith("emb_")]
    embedding_map = {
        str(row["protein_id"]): row[embedding_columns].to_numpy(dtype=np.float32)
        for _, row in embeddings.iterrows()
    }
    splits = pd.read_csv(SPLIT_PATH)
    return benchmark, targets, embedding_map, splits


def _stable_seed(text: str, seed: int) -> int:
    return seed + sum((index + 1) * ord(char) for index, char in enumerate(text)) % 1_000_003


def sample_rank_pairs(
    benchmark: pd.DataFrame,
    proteins: list[str],
    *,
    pairs_per_protein: int,
    seed: int,
    tie_tolerance: float = 1e-10,
) -> pd.DataFrame:
    """Sample deterministic within-protein pairs without materializing O(N^2) pairs."""
    output: list[dict[str, Any]] = []
    bins = [
        ("easy", 0.50, 0.40),
        ("medium", 0.20, 0.35),
        ("hard", 0.0, 0.25),
    ]
    for protein in sorted(proteins):
        group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
        scores = group["experimental_score"].to_numpy(dtype=float)
        ranks = pd.Series(scores).rank(method="average", pct=True).to_numpy()
        rng = np.random.default_rng(_stable_seed(protein, seed))
        selected: list[dict[str, Any]] = []
        target_counts = [int(round(pairs_per_protein * fraction)) for _, _, fraction in bins]
        target_counts[-1] = pairs_per_protein - sum(target_counts[:-1])
        for (label, lower, _), target_count in zip(bins, target_counts):
            attempts = 0
            while sum(row["difficulty"] == label for row in selected) < target_count and attempts < pairs_per_protein * 200:
                attempts += 1
                left = int(rng.integers(0, len(group)))
                right = int(rng.integers(0, len(group)))
                if left == right or abs(scores[left] - scores[right]) <= tie_tolerance:
                    continue
                rank_difference = abs(ranks[left] - ranks[right])
                if (label == "easy" and rank_difference < 0.50) or (
                    label == "medium" and not (0.20 <= rank_difference < 0.50)
                ) or (label == "hard" and not (0.0 < rank_difference < 0.20)):
                    continue
                selected.append(
                    {
                        "protein_id": protein,
                        "left_index": left,
                        "right_index": right,
                        "left_score": scores[left],
                        "right_score": scores[right],
                        "difficulty": label,
                    }
                )
            if sum(row["difficulty"] == label for row in selected) < target_count:
                raise RuntimeError(f"Unable to sample enough {label} pairs for {protein}")
        output.extend(selected[:pairs_per_protein])
    return pd.DataFrame(output)


def pairwise_ranking_loss(
    left_prediction: torch.Tensor,
    right_prediction: torch.Tensor,
    left_truth: torch.Tensor,
    right_truth: torch.Tensor,
    tie_tolerance: float = 1e-10,
) -> torch.Tensor:
    truth_delta = left_truth - right_truth
    valid = truth_delta.abs() > tie_tolerance
    if not torch.any(valid):
        return (left_prediction - right_prediction).sum() * 0.0
    sign = torch.sign(truth_delta[valid])
    prediction_delta = left_prediction[valid] - right_prediction[valid]
    return torch.nn.functional.softplus(-prediction_delta * sign).mean()


def _model_forward(
    model: nn.Module,
    model_name: str,
    embedding: torch.Tensor,
    target_features: torch.Tensor,
    candidate_features: torch.Tensor,
) -> torch.Tensor:
    if model_name == CandidateDNAOnly.model_name:
        return model(candidate_features)
    protein = embedding.reshape(1, -1).expand(candidate_features.shape[0], -1)
    target = target_features.unsqueeze(0).expand(candidate_features.shape[0], -1, -1)
    if model_name == ProteinCandidate.model_name:
        return model(protein, candidate_features)
    if model_name == TargetCandidateOnly.model_name:
        return model(target, candidate_features)
    if model_name == ProteinTargetCandidate.model_name:
        return model(protein, target, candidate_features)
    raise ValueError(f"Unknown model name: {model_name}")


def train_model(
    model: nn.Module,
    model_name: str,
    cache: dict[str, dict[str, Any]],
    pairs: pd.DataFrame,
    train_proteins: list[str],
    config: V05Config,
) -> tuple[nn.Module, pd.DataFrame, float]:
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay)
    history: list[dict[str, float]] = []
    start_time = time.perf_counter()
    for epoch in range(1, config.epochs + 1):
        model.train()
        losses = []
        ordered_pairs = pairs.sample(frac=1.0, random_state=config.seed + epoch).reset_index(drop=True)
        for protein in sorted(train_proteins):
            group_pairs = ordered_pairs.loc[ordered_pairs["protein_id"].eq(protein)]
            entry = cache[protein]
            for start in range(0, len(group_pairs), config.batch_size):
                batch = group_pairs.iloc[start : start + config.batch_size]
                left = torch.tensor(batch["left_index"].to_numpy(), dtype=torch.long)
                right = torch.tensor(batch["right_index"].to_numpy(), dtype=torch.long)
                candidate_features = entry["candidate_features"]
                left_prediction = _model_forward(
                    model,
                    model_name,
                    entry["embedding"],
                    entry["target_features"],
                    candidate_features[left],
                )
                right_prediction = _model_forward(
                    model,
                    model_name,
                    entry["embedding"],
                    entry["target_features"],
                    candidate_features[right],
                )
                left_truth = torch.tensor(batch["left_score"].to_numpy(), dtype=torch.float32)
                right_truth = torch.tensor(batch["right_score"].to_numpy(), dtype=torch.float32)
                loss = pairwise_ranking_loss(
                    left_prediction,
                    right_prediction,
                    left_truth,
                    right_truth,
                    config.tie_tolerance,
                )
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                losses.append(float(loss.detach()))
        history.append(
            {
                "model": model_name,
                "epoch": epoch,
                "mean_pairwise_loss": float(np.mean(losses)),
            }
        )
    return model, pd.DataFrame(history), time.perf_counter() - start_time


def build_cache(
    benchmark: pd.DataFrame,
    targets: pd.DataFrame,
    embedding_map: dict[str, np.ndarray],
    proteins: list[str],
) -> dict[str, dict[str, Any]]:
    target_map = targets.set_index("dbp_id")["primary_target"].to_dict()
    cache = {}
    for protein in sorted(proteins):
        group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
        candidates = group["candidate_dna"].astype(str).tolist()
        cache[protein] = {
            "candidate_dna": candidates,
            "truth": group["experimental_score"].to_numpy(dtype=float),
            "candidate_features": batch_rc_symmetric_one_hot(candidates),
            "target_features": target_window_features(target_map[protein]),
            "target": target_map[protein],
            "embedding": torch.from_numpy(embedding_map[protein]).float(),
        }
    return cache


def predict_model(
    model: nn.Module,
    model_name: str,
    entry: dict[str, Any],
    chunk_size: int = 512,
) -> np.ndarray:
    predictions = []
    candidates = entry["candidate_dna"]
    for start in range(0, len(candidates), chunk_size):
        end = start + chunk_size
        predictions.append(
            score_rc(
                model,
                model_name,
                candidates[start:end],
                protein_embedding=entry["embedding"],
                target=entry["target"],
            )
        )
    return np.concatenate(predictions)


def evaluate_predictions(
    benchmark: pd.DataFrame,
    predictions: dict[str, dict[str, np.ndarray]],
    test_proteins: list[str],
) -> pd.DataFrame:
    rows = []
    for protein in sorted(test_proteins):
        group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
        for model_name, per_protein in predictions.items():
            scored = group[["candidate_dna", "experimental_score"]].copy()
            scored["prediction"] = per_protein[protein]
            metrics = compute_ranking_metrics(scored, "experimental_score", "prediction")
            rows.append(
                {
                    "protein_id": protein,
                    "model": model_name,
                    "spearman": metrics.spearman,
                    "ndcg_1pct": metrics.ndcg_1pct,
                    "ndcg_5pct": metrics.ndcg_5pct,
                    "pairwise_accuracy": metrics.pairwise_accuracy,
                    "top1pct_recovery": metrics.top1pct_recovery,
                    "n_rc_classes": metrics.n_rc_classes,
                }
            )
    return pd.DataFrame(rows)


def evaluate_target_controls(
    benchmark: pd.DataFrame,
    targets: pd.DataFrame,
    test_proteins: list[str],
) -> pd.DataFrame:
    target_map = targets.set_index("dbp_id")["primary_target"].to_dict()
    rows = []
    for protein in sorted(test_proteins):
        group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
        for label, scorer in [
            ("TargetHamming", target_hamming_control),
            ("TargetEdit", target_edit_control),
            ("TargetKmerOverlap", target_kmer_overlap_control),
        ]:
            values = np.array([scorer(target_map[protein], candidate) for candidate in group["candidate_dna"]])
            scored = group[["candidate_dna", "experimental_score"]].copy()
            scored["prediction"] = values
            metrics = compute_ranking_metrics(scored, "experimental_score", "prediction")
            rows.append(
                {
                    "protein_id": protein,
                    "control": label,
                    "spearman": metrics.spearman,
                    "ndcg_1pct": metrics.ndcg_1pct,
                    "ndcg_5pct": metrics.ndcg_5pct,
                    "pairwise_accuracy": metrics.pairwise_accuracy,
                    "top1pct_recovery": metrics.top1pct_recovery,
                    "n_rc_classes": metrics.n_rc_classes,
                }
            )
    return pd.DataFrame(rows)


def run_smoke(config: V05Config) -> dict[str, pd.DataFrame]:
    set_seed(config.seed)
    benchmark, targets, embedding_map, splits = load_v05_data()
    fold = splits.loc[
        splits["split_name"].eq(config.smoke_split_name) & splits["fold_id"].eq(config.smoke_fold_id)
    ]
    if fold.empty:
        raise ValueError(f"Smoke fold not found: {config.smoke_split_name}/{config.smoke_fold_id}")
    train_proteins = sorted(fold.loc[fold["partition"].eq("train"), "dbp_id"].unique())
    test_proteins = sorted(fold.loc[fold["partition"].eq("test"), "dbp_id"].unique())
    if set(train_proteins) & set(test_proteins):
        raise AssertionError("Smoke split leaks protein IDs")
    pairs = sample_rank_pairs(
        benchmark,
        train_proteins,
        pairs_per_protein=config.pair_count_per_protein,
        seed=config.seed,
        tie_tolerance=config.tie_tolerance,
    )
    cache = build_cache(benchmark, targets, embedding_map, sorted(set(train_proteins + test_proteins)))
    protein_dim = len(next(iter(embedding_map.values())))
    model_factories = [
        (CandidateDNAOnly.model_name, lambda: CandidateDNAOnly(config.hidden_dim)),
        (ProteinCandidate.model_name, lambda: ProteinCandidate(protein_dim, config.hidden_dim)),
        (TargetCandidateOnly.model_name, lambda: TargetCandidateOnly(config.hidden_dim)),
        (ProteinTargetCandidate.model_name, lambda: ProteinTargetCandidate(protein_dim, config.hidden_dim)),
        (
            "M1c_ProteinCandidateCapacityMatched",
            lambda: ProteinCandidate(protein_dim, config.capacity_matched_hidden_dim),
        ),
    ]
    predictions: dict[str, dict[str, np.ndarray]] = {}
    history_parts = []
    runtime_rows = []
    parameter_rows = []
    health_rows = []
    for model_name, factory in model_factories:
        model = factory()
        trainable, frozen = model_parameter_counts(model)
        parameter_rows.append(
            {
                "model": model_name,
                "inputs": getattr(model, "inputs", "P,D"),
                "trainable_params": trainable,
                "frozen_params": frozen,
                "protein_representation": "none" if model_name == CandidateDNAOnly.model_name else "frozen ESM-2 t12 35M, 480 dimensions",
                "notes": "M1c is capacity-matched control; external frozen embeddings are not model parameters.",
            }
        )
        canonical_name = ProteinCandidate.model_name if model_name.startswith("M1c_") else model_name
        model, history, runtime = train_model(model, canonical_name, cache, pairs, train_proteins, config)
        history["reported_model"] = model_name
        history_parts.append(history)
        runtime_rows.append(
            {
                "model": model_name,
                "runtime_seconds": runtime,
                "training_proteins": "|".join(train_proteins),
                "test_proteins": "|".join(test_proteins),
            }
        )
        predictions[model_name] = {
            protein: predict_model(model, canonical_name, cache[protein]) for protein in test_proteins
        }
        all_predictions = np.concatenate(list(predictions[model_name].values()))
        last_loss = float(history["mean_pairwise_loss"].iloc[-1])
        health_rows.append(
            {
                "model": model_name,
                "prediction_variance": float(np.var(all_predictions)),
                "prediction_min": float(np.min(all_predictions)),
                "prediction_max": float(np.max(all_predictions)),
                "nan_inf_count": int(np.sum(~np.isfinite(all_predictions))),
                "first_epoch_loss": float(history["mean_pairwise_loss"].iloc[0]),
                "last_epoch_loss": last_loss,
                "loss_delta": last_loss - float(history["mean_pairwise_loss"].iloc[0]),
            }
        )
    evaluation = evaluate_predictions(benchmark, predictions, test_proteins)
    controls = evaluate_target_controls(benchmark, targets, test_proteins)
    smoke_info = pd.DataFrame(
        [
            {
                "split_name": config.smoke_split_name,
                "fold_id": config.smoke_fold_id,
                "train_proteins": "|".join(train_proteins),
                "test_proteins": "|".join(test_proteins),
                "n_train_proteins": len(train_proteins),
                "n_test_proteins": len(test_proteins),
                "n_rc_units_per_protein": 8192,
                "training_pair_count": len(pairs),
                "status": "NOT PRIMARY SCIENTIFIC RESULT",
            }
        ]
    )
    return {
        "evaluation": evaluation,
        "controls": controls,
        "pairs": pairs,
        "history": pd.concat(history_parts, ignore_index=True),
        "runtime": pd.DataFrame(runtime_rows),
        "parameters": pd.DataFrame(parameter_rows),
        "health": pd.DataFrame(health_rows),
        "smoke_info": smoke_info,
    }


def write_smoke_outputs(config: V05Config, outputs: dict[str, pd.DataFrame]) -> None:
    results = ensure_dir(ROOT / "results" / "v0_5")
    logs = ensure_dir(ROOT / "logs" / "v0_5")
    metadata = ensure_dir(ROOT / "metadata" / "v0_5")
    outputs["evaluation"].to_csv(results / "smoke_test_results.csv", index=False)
    outputs["controls"].to_csv(results / "smoke_target_control_results.csv", index=False)
    outputs["parameters"].to_csv(results / "model_parameter_counts.csv", index=False)
    outputs["history"].to_csv(results / "smoke_loss_history.csv", index=False)
    outputs["runtime"].to_csv(results / "smoke_runtime.csv", index=False)
    outputs["health"].to_csv(results / "smoke_training_health.csv", index=False)
    outputs["smoke_info"].to_csv(results / "smoke_split_info.csv", index=False)
    outputs["pairs"].to_csv(results / "smoke_training_pairs.csv", index=False)
    config_path = metadata / "v0_5_model_config.json"
    config_path.write_text(json.dumps(config.as_dict(), indent=2), encoding="utf-8")
    (logs / "01_model_smoke.log").write_text(
        outputs["smoke_info"].to_string(index=False)
        + "\n\n"
        + outputs["runtime"].to_string(index=False)
        + "\n\n"
        + outputs["evaluation"].to_string(index=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    config = V05Config()
    result = run_smoke(config)
    write_smoke_outputs(config, result)
    print(result["smoke_info"].to_string(index=False))
    print(result["evaluation"].to_string(index=False))
