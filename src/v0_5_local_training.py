from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from torch import nn

from src.utils import ensure_dir, project_root
from src.v0_4_evaluation import compute_ranking_metrics
from src.v0_5_local_models import (
    LocalProteinCandidate,
    LocalProteinCandidateCapacityMatched,
    LocalProteinTargetCandidate,
    batch_rc_oriented_dna,
    local_model_parameter_counts,
    score_local_rc,
    target_window_oriented_features,
)
from src.v0_5_training import pairwise_ranking_loss, sample_rank_pairs, set_seed


ROOT = project_root()
BENCHMARK_PATH = ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet"
TARGET_PATH = ROOT / "metadata" / "v0_5" / "designed_target_manifest_v0_5.csv"
RESIDUE_EMBEDDING_PATH = (
    ROOT / "data" / "interim" / "v0_5_local" / "designed_residue_embeddings_esm2_t12_35M_UR50D.parquet"
)
SPLIT_PATH = ROOT / "metadata" / "v0_5" / "v0_5_split_manifest.csv"
OLD_PRIMARY_PATH = ROOT / "results" / "v0_5" / "primary_per_protein_results.csv"
LOCAL_CONFIG_PATH = ROOT / "metadata" / "v0_5_local" / "local_model_config.json"


@dataclass(frozen=True)
class LocalConfig:
    seed: int = 42
    smoke_split_name: str = "protein_cluster_loco"
    smoke_fold_id: str = "protein_cluster_loco_fold_2"
    hidden_dim: int = 24
    capacity_matched_head_hidden_dim: int = 64
    learning_rate: float = 0.01
    weight_decay: float = 1e-4
    epochs: int = 18
    pair_count_per_protein: int = 512
    batch_size: int = 128
    tie_tolerance: float = 1e-10

    @classmethod
    def from_json(cls, path: Path = LOCAL_CONFIG_PATH) -> "LocalConfig":
        values = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            seed=int(values["seed"]),
            smoke_split_name=str(values["smoke_split_name"]),
            smoke_fold_id=str(values["smoke_fold_id"]),
            hidden_dim=int(values["hidden_dim"]),
            capacity_matched_head_hidden_dim=int(values["capacity_matched_head_hidden_dim"]),
            learning_rate=float(values["learning_rate"]),
            weight_decay=float(values["weight_decay"]),
            epochs=int(values["epochs"]),
            pair_count_per_protein=int(values["pair_count_per_protein"]),
            batch_size=int(values["batch_size"]),
            tie_tolerance=float(values["tie_tolerance"]),
        )


def load_local_data() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, torch.Tensor], pd.DataFrame]:
    benchmark = pd.read_parquet(BENCHMARK_PATH).rename(
        columns={
            "canonical_7mer": "candidate_dna",
            "experimental_escore_consensus": "experimental_score",
        }
    )
    benchmark = benchmark[["protein_id", "candidate_dna", "experimental_score"]].copy()
    targets = pd.read_csv(TARGET_PATH)
    residues = pd.read_parquet(RESIDUE_EMBEDDING_PATH)
    embedding_columns = [column for column in residues.columns if column.startswith("emb_")]
    residue_map: dict[str, torch.Tensor] = {}
    for protein_id, group in residues.groupby("protein_id", sort=False):
        group = group.sort_values("residue_index")
        residue_map[str(protein_id)] = torch.from_numpy(
            group[embedding_columns].to_numpy(dtype=np.float32, copy=True)
        )
    splits = pd.read_csv(SPLIT_PATH)
    return benchmark, targets, residue_map, splits


def build_local_cache(
    benchmark: pd.DataFrame,
    targets: pd.DataFrame,
    residue_map: dict[str, torch.Tensor],
    proteins: list[str],
) -> dict[str, dict[str, Any]]:
    target_map = targets.set_index("dbp_id")["primary_target"].to_dict()
    cache: dict[str, dict[str, Any]] = {}
    for protein in sorted(proteins):
        if protein not in residue_map:
            raise KeyError(f"Missing residue embedding for {protein}")
        group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
        candidates = group["candidate_dna"].astype(str).tolist()
        cache[protein] = {
            "candidate_dna": candidates,
            "truth": group["experimental_score"].to_numpy(dtype=float),
            "candidate_features": batch_rc_oriented_dna(candidates),
            "target": str(target_map[protein]),
            "target_features": target_window_oriented_features(str(target_map[protein])),
            "residue_embeddings": residue_map[protein].float(),
        }
    return cache


def _forward(
    model: nn.Module,
    model_name: str,
    entry: dict[str, Any],
    candidate_features: torch.Tensor,
) -> torch.Tensor:
    residues = entry["residue_embeddings"]
    if model_name in {
        LocalProteinCandidate.model_name,
        LocalProteinCandidateCapacityMatched.model_name,
    }:
        return model(residues, candidate_features)
    if model_name == LocalProteinTargetCandidate.model_name:
        return model(residues, entry["target_features"], candidate_features)
    raise ValueError(f"Unknown local model name: {model_name}")


def train_local_model(
    model: nn.Module,
    model_name: str,
    cache: dict[str, dict[str, Any]],
    pairs: pd.DataFrame,
    train_proteins: list[str],
    config: LocalConfig,
) -> tuple[nn.Module, pd.DataFrame, float]:
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay)
    history: list[dict[str, float]] = []
    start_time = time.perf_counter()
    for epoch in range(1, config.epochs + 1):
        model.train()
        losses: list[float] = []
        ordered_pairs = pairs.sample(frac=1.0, random_state=config.seed + epoch).reset_index(drop=True)
        for protein in sorted(train_proteins):
            protein_pairs = ordered_pairs.loc[ordered_pairs["protein_id"].eq(protein)]
            entry = cache[protein]
            for start in range(0, len(protein_pairs), config.batch_size):
                batch = protein_pairs.iloc[start : start + config.batch_size]
                left = torch.tensor(batch["left_index"].to_numpy(), dtype=torch.long)
                right = torch.tensor(batch["right_index"].to_numpy(), dtype=torch.long)
                left_prediction = _forward(model, model_name, entry, entry["candidate_features"][left])
                right_prediction = _forward(model, model_name, entry, entry["candidate_features"][right])
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


def predict_local_model(
    model: nn.Module,
    model_name: str,
    entry: dict[str, Any],
    chunk_size: int = 256,
) -> np.ndarray:
    predictions = []
    candidates = entry["candidate_dna"]
    for start in range(0, len(candidates), chunk_size):
        predictions.append(
            score_local_rc(
                model,
                model_name,
                entry["residue_embeddings"],
                candidates[start : start + chunk_size],
                target=entry["target"],
            )
        )
    return np.concatenate(predictions)


def _prediction_correlation(first: np.ndarray, second: np.ndarray) -> float:
    if np.std(first) == 0.0 or np.std(second) == 0.0:
        return np.nan
    return float(np.corrcoef(first, second)[0, 1])


def _spearman(truth: np.ndarray, prediction: np.ndarray) -> float:
    if np.unique(truth).size < 2 or np.unique(prediction).size < 2:
        return np.nan
    return float(spearmanr(truth, prediction).statistic)


def _attention_stats(attention: torch.Tensor) -> dict[str, float]:
    attention = attention.detach().cpu().numpy()
    epsilon = 1e-12
    length = attention.shape[-1]
    entropy = -np.sum(attention * np.log(attention + epsilon), axis=-1) / np.log(max(length, 2))
    variation = float(np.std(attention, axis=0).mean())
    top_residues = np.argmax(attention, axis=-1).reshape(-1)
    counts = np.bincount(top_residues, minlength=length)
    return {
        "attention_entropy_normalized": float(np.mean(entropy)),
        "candidate_attention_variation_mean_std": variation,
        "top_residue_max_fraction": float(np.max(counts) / len(top_residues)),
        "top_residue_unique_count": int(np.count_nonzero(counts)),
    }


def attention_diagnostics(
    models: dict[str, nn.Module],
    cache: dict[str, dict[str, Any]],
    test_proteins: list[str],
    sample_size: int = 256,
) -> pd.DataFrame:
    rows = []
    for model_name, model in models.items():
        for protein in sorted(test_proteins):
            entry = cache[protein]
            candidates = entry["candidate_features"][:sample_size]
            model.eval()
            with torch.no_grad():
                if model_name == LocalProteinTargetCandidate.model_name:
                    _, attention, _ = model.forward_with_attention(
                        entry["residue_embeddings"],
                        entry["target_features"],
                        candidates,
                    )
                else:
                    _, attention = model.forward_with_attention(entry["residue_embeddings"], candidates)
            stats = _attention_stats(attention)
            rows.append(
                {
                    "model": model_name,
                    "dbp_id": protein,
                    "n_candidates_sampled": min(sample_size, len(entry["candidate_dna"])),
                    **stats,
                    "interpretation": "mechanism sanity diagnostic; not a biological contact map",
                }
            )
    return pd.DataFrame(rows)


def shuffle_diagnostics(
    models: dict[str, nn.Module],
    cache: dict[str, dict[str, Any]],
    all_proteins: list[str],
    test_proteins: list[str],
) -> pd.DataFrame:
    ordered = sorted(all_proteins)
    next_protein = {protein: ordered[(index + 1) % len(ordered)] for index, protein in enumerate(ordered)}
    rows = []
    for model_name, model in models.items():
        for protein in sorted(test_proteins):
            entry = cache[protein]
            original = predict_local_model(model, model_name, entry)
            original_spearman = _spearman(entry["truth"], original)
            shuffled_p_entry = dict(entry)
            shuffled_p_entry["residue_embeddings"] = cache[next_protein[protein]]["residue_embeddings"]
            shuffled_p = predict_local_model(model, model_name, shuffled_p_entry)
            rows.append(
                {
                    "model": model_name,
                    "dbp_id": protein,
                    "shuffle_type": "protein",
                    "permutation_source": next_protein[protein],
                    "prediction_correlation": _prediction_correlation(original, shuffled_p),
                    "original_spearman": original_spearman,
                    "shuffled_spearman": _spearman(entry["truth"], shuffled_p),
                    "spearman_delta": _spearman(entry["truth"], shuffled_p) - original_spearman,
                    "mean_abs_score_change": float(np.mean(np.abs(original - shuffled_p))),
                    "retrained": False,
                }
            )
            if model_name == LocalProteinTargetCandidate.model_name:
                shuffled_target_entry = dict(entry)
                shuffled_target_entry["target"] = cache[next_protein[protein]]["target"]
                shuffled_target_entry["target_features"] = cache[next_protein[protein]]["target_features"]
                shuffled_t = predict_local_model(model, model_name, shuffled_target_entry)
                rows.append(
                    {
                        "model": model_name,
                        "dbp_id": protein,
                        "shuffle_type": "target",
                        "permutation_source": next_protein[protein],
                        "prediction_correlation": _prediction_correlation(original, shuffled_t),
                        "original_spearman": original_spearman,
                        "shuffled_spearman": _spearman(entry["truth"], shuffled_t),
                        "spearman_delta": _spearman(entry["truth"], shuffled_t) - original_spearman,
                        "mean_abs_score_change": float(np.mean(np.abs(original - shuffled_t))),
                        "retrained": False,
                    }
                )
    return pd.DataFrame(rows)


def run_local_smoke(config: LocalConfig | None = None) -> dict[str, pd.DataFrame]:
    config = LocalConfig() if config is None else config
    set_seed(config.seed)
    torch.set_num_threads(2)
    benchmark, targets, residue_map, splits = load_local_data()
    fold = splits.loc[
        splits["split_name"].eq(config.smoke_split_name) & splits["fold_id"].eq(config.smoke_fold_id)
    ]
    if fold.empty:
        raise ValueError(f"Smoke fold not found: {config.smoke_split_name}/{config.smoke_fold_id}")
    train_proteins = sorted(fold.loc[fold["partition"].eq("train"), "dbp_id"].unique())
    test_proteins = sorted(fold.loc[fold["partition"].eq("test"), "dbp_id"].unique())
    if set(train_proteins) & set(test_proteins):
        raise AssertionError("Local smoke split leaks protein IDs")
    pairs = sample_rank_pairs(
        benchmark,
        train_proteins,
        pairs_per_protein=config.pair_count_per_protein,
        seed=config.seed,
        tie_tolerance=config.tie_tolerance,
    )
    all_proteins = sorted(set(train_proteins + test_proteins))
    cache = build_local_cache(benchmark, targets, residue_map, all_proteins)
    model_factories = [
        (
            LocalProteinCandidate.model_name,
            lambda: LocalProteinCandidate(protein_dim=480, hidden_dim=config.hidden_dim),
        ),
        (
            LocalProteinCandidateCapacityMatched.model_name,
            lambda: LocalProteinCandidateCapacityMatched(
                protein_dim=480,
                hidden_dim=config.hidden_dim,
                head_hidden_dim=config.capacity_matched_head_hidden_dim,
            ),
        ),
        (
            LocalProteinTargetCandidate.model_name,
            lambda: LocalProteinTargetCandidate(protein_dim=480, hidden_dim=config.hidden_dim),
        ),
    ]
    trained_models: dict[str, nn.Module] = {}
    evaluation_rows = []
    health_rows = []
    parameter_rows = []
    history_rows = []
    runtime_rows = []
    for model_name, factory in model_factories:
        model = factory()
        trainable, frozen = local_model_parameter_counts(model)
        parameter_rows.append(
            {
                "model": model_name,
                "inputs": getattr(model, "inputs", "P,D"),
                "trainable_params": trainable,
                "frozen_params": frozen,
                "frozen_protein_embedding_dim": 480,
                "notes": "ESM residue embeddings are external frozen inputs; not trainable model parameters.",
            }
        )
        model, history, runtime = train_local_model(model, model_name, cache, pairs, train_proteins, config)
        trained_models[model_name] = model
        history_rows.append(history)
        runtime_rows.append(
            {
                "model": model_name,
                "runtime_seconds": runtime,
                "training_proteins": "|".join(train_proteins),
                "test_proteins": "|".join(test_proteins),
            }
        )
        train_predictions = {
            protein: predict_local_model(model, model_name, cache[protein]) for protein in train_proteins
        }
        test_predictions = {
            protein: predict_local_model(model, model_name, cache[protein]) for protein in test_proteins
        }
        all_prediction_values = np.concatenate(list(test_predictions.values()))
        for protein in all_proteins:
            predictions = train_predictions.get(protein, test_predictions.get(protein))
            group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
            metrics = compute_ranking_metrics(
                pd.DataFrame({"experimental_score": group["experimental_score"], "prediction": predictions}),
                "experimental_score",
                "prediction",
            )
            evaluation_rows.append(
                {
                    "model": model_name,
                    "dbp_id": protein,
                    "partition": "train" if protein in train_proteins else "test",
                    "spearman": metrics.spearman,
                    "ndcg_1pct": metrics.ndcg_1pct,
                    "ndcg_5pct": metrics.ndcg_5pct,
                    "pairwise_accuracy": metrics.pairwise_accuracy,
                    "top1pct_recovery": metrics.top1pct_recovery,
                    "n_rc_units": metrics.n_rc_classes,
                }
            )
        last_loss = float(history["mean_pairwise_loss"].iloc[-1])
        first_loss = float(history["mean_pairwise_loss"].iloc[0])
        train_values = [
            row["spearman"]
            for row in evaluation_rows
            if row["model"] == model_name and row["partition"] == "train" and np.isfinite(row["spearman"])
        ]
        train_spearman = float(np.median(train_values)) if train_values else np.nan
        health_rows.append(
            {
                "model": model_name,
                "training_pair_count": len(pairs),
                "first_epoch_loss": first_loss,
                "final_epoch_loss": last_loss,
                "loss_delta": last_loss - first_loss,
                "train_macro_median_spearman": train_spearman,
                "test_prediction_variance": float(np.var(all_prediction_values)),
                "test_prediction_min": float(np.min(all_prediction_values)),
                "test_prediction_max": float(np.max(all_prediction_values)),
                "nan_inf_count": int(np.sum(~np.isfinite(all_prediction_values))),
                "runtime_seconds": runtime,
                "status": "complete",
            }
        )
    evaluation = pd.DataFrame(evaluation_rows)
    old = pd.read_csv(OLD_PRIMARY_PATH).set_index("dbp_id")
    wide_rows = []
    for protein in test_proteins:
        local = evaluation.loc[(evaluation["dbp_id"].eq(protein)) & evaluation["partition"].eq("test")]
        row = {
            "dbp_id": protein,
            "old_M1c": float(old.loc[protein, "M1c"]),
            "old_M2": float(old.loc[protein, "M2"]),
            "old_M3": float(old.loc[protein, "M3"]),
            "L1": float(local.loc[local["model"].eq(LocalProteinCandidate.model_name), "spearman"].iloc[0]),
            "L1c": float(
                local.loc[local["model"].eq(LocalProteinCandidateCapacityMatched.model_name), "spearman"].iloc[0]
            ),
            "L2": float(local.loc[local["model"].eq(LocalProteinTargetCandidate.model_name), "spearman"].iloc[0]),
            "context_note": "old columns are frozen v0.5 primary seed-mean values; local columns are Phase 5A smoke only",
        }
        wide_rows.append(row)
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
                "exposure_status": "development-exposed-for-local-model",
            }
        ]
    )
    return {
        "smoke_results": pd.DataFrame(wide_rows),
        "evaluation": evaluation,
        "parameters": pd.DataFrame(parameter_rows),
        "health": pd.DataFrame(health_rows),
        "history": pd.concat(history_rows, ignore_index=True),
        "runtime": pd.DataFrame(runtime_rows),
        "pairs": pairs,
        "smoke_info": smoke_info,
        "shuffle": shuffle_diagnostics(trained_models, cache, all_proteins, test_proteins),
        "attention": attention_diagnostics(trained_models, cache, test_proteins),
    }


def write_local_smoke_outputs(outputs: dict[str, pd.DataFrame]) -> None:
    results = ensure_dir(ROOT / "results" / "v0_5_local")
    logs = ensure_dir(ROOT / "logs" / "v0_5_local")
    outputs["smoke_results"].to_csv(results / "smoke_results.csv", index=False)
    outputs["evaluation"].to_csv(results / "smoke_per_protein_results.csv", index=False)
    outputs["parameters"].to_csv(results / "model_parameter_counts.csv", index=False)
    outputs["health"].to_csv(results / "smoke_training_health.csv", index=False)
    outputs["history"].to_csv(results / "smoke_loss_history.csv", index=False)
    outputs["runtime"].to_csv(results / "smoke_runtime.csv", index=False)
    outputs["pairs"].to_csv(results / "smoke_training_pairs.csv", index=False)
    outputs["smoke_info"].to_csv(results / "smoke_split_info.csv", index=False)
    outputs["shuffle"].to_csv(results / "shuffle_diagnostics.csv", index=False)
    outputs["attention"].to_csv(results / "attention_diagnostics.csv", index=False)
    log = [
        "Phase 5A local interaction smoke",
        outputs["smoke_info"].to_string(index=False),
        outputs["parameters"].to_string(index=False),
        outputs["health"].to_string(index=False),
        outputs["smoke_results"].to_string(index=False),
        outputs["shuffle"].to_string(index=False),
        outputs["attention"].to_string(index=False),
    ]
    (logs / "01_local_model_smoke.log").write_text("\n\n".join(log), encoding="utf-8")


if __name__ == "__main__":
    result = run_local_smoke(LocalConfig.from_json())
    write_local_smoke_outputs(result)
    print(result["smoke_results"].to_string(index=False))
