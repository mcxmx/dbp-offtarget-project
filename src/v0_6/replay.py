from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd
import torch

from src.utils import project_root
from src.v0_4_evaluation import compute_ranking_metrics
from src.v0_5_primary_evaluation import get_fold_partitions
from src.v0_5_training import load_v05_data, sample_rank_pairs, set_seed

from .models import (
    CandidateDNAOnlyV06,
    ProteinCandidateV06,
    ProteinTargetCandidateV06,
    TargetCandidateOnlyV06,
    batch_canonical_one_hot,
    model_parameter_counts_v06,
    score_rc_v06,
    target_window_features_v06,
)
from .audit import DESIGNED_IDS, MODEL_ORDER


ROOT = project_root()
RESULTS = ROOT / "results" / "v0_6_representation"
SEEDS = (17, 29, 43)
V05_PRIMARY = ROOT / "results" / "v0_5"


@dataclass(frozen=True)
class V06ReplayConfig:
    seed: int = 42
    hidden_dim: int = 32
    capacity_matched_hidden_dim: int = 40
    learning_rate: float = 0.01
    weight_decay: float = 1e-4
    epochs: int = 18
    pair_count_per_protein: int = 512
    batch_size: int = 128
    tie_tolerance: float = 1e-10


def build_cache_v06(
    benchmark: pd.DataFrame,
    targets: pd.DataFrame,
    embedding_map: dict[str, np.ndarray],
    proteins: Iterable[str],
) -> dict[str, dict[str, Any]]:
    target_map = targets.set_index("dbp_id")["primary_target"].to_dict()
    cache: dict[str, dict[str, Any]] = {}
    for protein in sorted(proteins):
        group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
        candidates = group["candidate_dna"].astype(str).tolist()
        cache[protein] = {
            "candidate_dna": candidates,
            "truth": group["experimental_score"].to_numpy(dtype=float),
            "candidate_features": batch_canonical_one_hot(candidates),
            "target_features": target_window_features_v06(target_map[protein]),
            "target": target_map[protein],
            "embedding": torch.from_numpy(np.array(embedding_map[protein], dtype=np.float32, copy=True)).float(),
        }
    return cache


def model_factory_v06(label: str, protein_dim: int, config: V06ReplayConfig):
    if label == "M0":
        return CandidateDNAOnlyV06(config.hidden_dim), CandidateDNAOnlyV06.model_name
    if label == "M1":
        return ProteinCandidateV06(protein_dim, config.hidden_dim), ProteinCandidateV06.model_name
    if label == "M1c":
        return ProteinCandidateV06(protein_dim, config.capacity_matched_hidden_dim), ProteinCandidateV06.model_name
    if label == "M2":
        return TargetCandidateOnlyV06(config.hidden_dim), TargetCandidateOnlyV06.model_name
    if label == "M3":
        return ProteinTargetCandidateV06(protein_dim, config.hidden_dim), ProteinTargetCandidateV06.model_name
    raise ValueError(label)


def _forward(model, name: str, entry: dict[str, Any], indices: np.ndarray) -> torch.Tensor:
    candidate = entry["candidate_features"][indices]
    n = len(indices)
    if name == CandidateDNAOnlyV06.model_name:
        return model(candidate)
    protein = entry["embedding"].reshape(1, -1).expand(n, -1)
    if name == ProteinCandidateV06.model_name:
        return model(protein, candidate)
    target = entry["target_features"].unsqueeze(0).expand(n, -1, -1)
    if name == TargetCandidateOnlyV06.model_name:
        return model(target, candidate)
    return model(protein, target, candidate)


def train_model_v06(
    model: torch.nn.Module,
    model_name: str,
    cache: dict[str, dict[str, Any]],
    pairs: pd.DataFrame,
    train_proteins: list[str],
    config: V06ReplayConfig,
) -> tuple[torch.nn.Module, pd.DataFrame]:
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay)
    history = []
    for epoch in range(1, config.epochs + 1):
        losses = []
        ordered = pairs.sample(frac=1.0, random_state=config.seed + epoch).reset_index(drop=True)
        for protein in sorted(train_proteins):
            entry = cache[protein]
            group = ordered.loc[ordered["protein_id"].eq(protein)]
            for start in range(0, len(group), config.batch_size):
                batch = group.iloc[start : start + config.batch_size]
                left = batch["left_index"].to_numpy(dtype=int)
                right = batch["right_index"].to_numpy(dtype=int)
                left_prediction = _forward(model, model_name, entry, left)
                right_prediction = _forward(model, model_name, entry, right)
                truth_delta = torch.tensor(batch["left_score"].to_numpy() - batch["right_score"].to_numpy(), dtype=torch.float32)
                valid = truth_delta.abs() > config.tie_tolerance
                sign = torch.sign(truth_delta[valid])
                loss = torch.nn.functional.softplus(
                    -(left_prediction[valid] - right_prediction[valid]) * sign
                ).mean()
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                losses.append(float(loss.detach()))
        history.append({"model": model_name, "epoch": epoch, "mean_pairwise_loss": float(np.mean(losses))})
    return model, pd.DataFrame(history)


def predict_v06(model, name: str, entry: dict[str, Any], chunk_size: int = 512) -> np.ndarray:
    output = []
    for start in range(0, len(entry["candidate_dna"]), chunk_size):
        end = start + chunk_size
        output.append(
            score_rc_v06(
                model,
                name,
                entry["candidate_dna"][start:end],
                protein_embedding=entry["embedding"],
                target=entry["target"],
            )
        )
    return np.concatenate(output)


def _metrics(truth: np.ndarray, prediction: np.ndarray) -> dict[str, Any]:
    frame = pd.DataFrame({"experimental_score": truth, "prediction": prediction})
    result = compute_ranking_metrics(frame, "experimental_score", "prediction")
    return {
        "spearman": result.spearman,
        "ndcg_1pct": result.ndcg_1pct,
        "ndcg_5pct": result.ndcg_5pct,
        "pairwise_accuracy": result.pairwise_accuracy,
        "top1pct_recovery": result.top1pct_recovery,
        "n_rc_classes": result.n_rc_classes,
    }


def deterministic_permutation(ids: Iterable[str]) -> dict[str, str]:
    ordered = sorted(set(ids))
    return {protein: ordered[(index + 1) % len(ordered)] for index, protein in enumerate(ordered)}


def run_v06_replay(config: V06ReplayConfig | None = None) -> dict[str, pd.DataFrame]:
    config = config or V06ReplayConfig()
    benchmark, targets, embedding_map, splits = load_v05_data()
    protein_dim = len(next(iter(embedding_map.values())))
    target_map = targets.set_index("dbp_id")["primary_target"].to_dict()
    protein_perm = deterministic_permutation(DESIGNED_IDS)
    target_perm = deterministic_permutation(DESIGNED_IDS)
    result_rows: list[dict[str, Any]] = []
    health_rows: list[dict[str, Any]] = []
    protein_shuffle_rows: list[dict[str, Any]] = []
    target_shuffle_rows: list[dict[str, Any]] = []
    parameter_rows: list[dict[str, Any]] = []
    for fold_id in sorted(splits.loc[splits["split_name"].eq("protein_cluster_loco"), "fold_id"].unique()):
        train_proteins, test_proteins = get_fold_partitions(splits, "protein_cluster_loco", fold_id)
        cache = build_cache_v06(benchmark, targets, embedding_map, train_proteins + test_proteins)
        for seed in SEEDS:
            set_seed(seed)
            pairs = sample_rank_pairs(benchmark, train_proteins, pairs_per_protein=512, seed=seed)
            for label in MODEL_ORDER:
                set_seed(seed)
                model, name = model_factory_v06(label, protein_dim, config)
                parameter_rows.append(
                    {
                        "model": label,
                        "trainable_params": model_parameter_counts_v06(model)[0],
                        "representation": "v0.6 canonical RC one-hot + ordered target windows",
                    }
                )
                model, history = train_model_v06(
                    model,
                    name,
                    cache,
                    pairs,
                    train_proteins,
                    V06ReplayConfig(**{**config.__dict__, "seed": int(seed)}),
                )
                test_predictions = {}
                for protein in test_proteins:
                    prediction = predict_v06(model, name, cache[protein])
                    test_predictions[protein] = prediction
                    metric = _metrics(cache[protein]["truth"], prediction)
                    result_rows.append(
                        {
                            "split_type": "protein_cluster_loco",
                            "fold_id": fold_id,
                            "seed": seed,
                            "dbp_id": protein,
                            "model": label,
                            "spearman": metric["spearman"],
                            "ndcg_1pct": metric["ndcg_1pct"],
                            "ndcg_5pct": metric["ndcg_5pct"],
                            "pairwise_accuracy": metric["pairwise_accuracy"],
                            "top1pct_recovery": metric["top1pct_recovery"],
                            "n_rc_classes": metric["n_rc_classes"],
                            "prediction_variance": float(np.var(prediction)),
                            "training_pair_count": len(pairs),
                            "development_exposed": True,
                            "status": "diagnostic_replay_all_designed_exposed",
                        }
                    )
                    if label == "M3":
                        original = prediction
                        shuffled_p = score_rc_v06(
                            model,
                            name,
                            cache[protein]["candidate_dna"],
                            protein_embedding=cache[protein_perm[protein]]["embedding"],
                            target=cache[protein]["target"],
                        )
                        shuffled_t = score_rc_v06(
                            model,
                            name,
                            cache[protein]["candidate_dna"],
                            protein_embedding=cache[protein]["embedding"],
                            target=target_map[target_perm[protein]],
                        )
                        protein_shuffle_rows.append(
                            {
                                "split_type": "protein_cluster_loco",
                                "fold_id": fold_id,
                                "seed": seed,
                                "dbp_id": protein,
                                "replacement_protein_id": protein_perm[protein],
                                "prediction_correlation": float(np.corrcoef(original, shuffled_p)[0, 1]),
                                "mean_abs_score_change": float(np.mean(np.abs(original - shuffled_p))),
                                "protein_conditioning_effect_size": float(1.0 - np.corrcoef(original, shuffled_p)[0, 1]),
                                "retrained": False,
                            }
                        )
                        target_shuffle_rows.append(
                            {
                                "split_type": "protein_cluster_loco",
                                "fold_id": fold_id,
                                "seed": seed,
                                "dbp_id": protein,
                                "replacement_target_dbp_id": target_perm[protein],
                                "prediction_correlation": float(np.corrcoef(original, shuffled_t)[0, 1]),
                                "mean_abs_score_change": float(np.mean(np.abs(original - shuffled_t))),
                                "target_conditioning_effect_size": float(1.0 - np.corrcoef(original, shuffled_t)[0, 1]),
                                "retrained": False,
                            }
                        )
                all_prediction = np.concatenate(list(test_predictions.values()))
                health_rows.append(
                    {
                        "split_type": "protein_cluster_loco",
                        "fold_id": fold_id,
                        "seed": seed,
                        "model": label,
                        "first_epoch_loss": history["mean_pairwise_loss"].iloc[0],
                        "last_epoch_loss": history["mean_pairwise_loss"].iloc[-1],
                        "prediction_variance": float(np.var(all_prediction)),
                        "nan_inf_count": int(np.sum(~np.isfinite(all_prediction))),
                    }
                )
    result = pd.DataFrame(result_rows)
    # Sequence-only k-mer baseline is copied as a declared historical baseline,
    # not recalculated from any v0.6 test outcome.
    kmer = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "designed_dbp_sequence_baseline_rc_aware.csv")
    kmer = kmer.loc[kmer["metric"].eq("kmer3_jaccard_to_paper_motif_rc_aware"), ["protein_id", "spearman"]].rename(columns={"protein_id": "dbp_id"})
    kmer["model"] = "sequence_kmer3_rc_aware_baseline"
    kmer["split_type"] = "historical_v0_3_1_no_training"
    kmer["seed"] = np.nan
    result = pd.concat([result, kmer[["split_type", "seed", "dbp_id", "model", "spearman"]]], ignore_index=True, sort=False)
    RESULTS.mkdir(parents=True, exist_ok=True)
    result.to_csv(RESULTS / "primary_results.csv", index=False)
    pd.DataFrame(health_rows).to_csv(RESULTS / "training_health.csv", index=False)
    pd.DataFrame(protein_shuffle_rows).to_csv(RESULTS / "protein_shuffle.csv", index=False)
    pd.DataFrame(target_shuffle_rows).to_csv(RESULTS / "target_shuffle.csv", index=False)
    pd.DataFrame(parameter_rows).drop_duplicates().to_csv(RESULTS / "model_parameters.csv", index=False)
    return {
        "primary": result,
        "health": pd.DataFrame(health_rows),
        "protein_shuffle": pd.DataFrame(protein_shuffle_rows),
        "target_shuffle": pd.DataFrame(target_shuffle_rows),
    }


if __name__ == "__main__":
    outputs = run_v06_replay()
    for frame in outputs.values():
        print(frame.head(20).to_string(index=False))
