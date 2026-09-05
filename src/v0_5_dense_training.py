from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.utils import ensure_dir, project_root
from src.v0_4_evaluation import compute_ranking_metrics
from src.v0_5_dense_sampling import audit_pair_sampling, sample_dense_pairs
from src.v0_5_models import CandidateDNAOnly, ProteinCandidate, ProteinTargetCandidate, TargetCandidateOnly
from src.v0_5_training import (
    V05Config,
    build_cache,
    evaluate_predictions,
    load_v05_data,
    predict_model,
    sample_rank_pairs,
    set_seed,
    train_model,
)


ROOT = project_root()
DENSE_CONFIG_PATH = ROOT / "metadata" / "v0_5_dense" / "dense_training_config.json"
OUTPUT_DIR = ROOT / "results" / "v0_5_dense"
LOG_DIR = ROOT / "logs" / "v0_5_dense"

MODEL_LABELS = {
    "M0": CandidateDNAOnly.model_name,
    "M1c": "M1c_ProteinCandidateCapacityMatched",
    "M2": TargetCandidateOnly.model_name,
    "M3": ProteinTargetCandidate.model_name,
}


def select_untouched_smoke_fold(
    splits: pd.DataFrame,
    *,
    still_untouched: set[str],
    split_name: str = "protein_cluster_loco",
) -> tuple[str, list[str], list[str]]:
    folds = splits.loc[splits["split_name"].eq(split_name), "fold_id"].drop_duplicates().tolist()
    for fold_id in folds:
        fold = splits.loc[splits["fold_id"].eq(fold_id)]
        test = sorted(fold.loc[fold["partition"].eq("test"), "dbp_id"].unique())
        train = sorted(fold.loc[fold["partition"].eq("train"), "dbp_id"].unique())
        if test and set(test).issubset(still_untouched):
            return fold_id, train, test
    raise ValueError("No fold consists entirely of still-untouched proteins")


def dense_config_from_registered() -> V05Config:
    values = __import__("json").loads(DENSE_CONFIG_PATH.read_text(encoding="utf-8"))
    return replace(
        V05Config(),
        seed=int(values["seed"]),
        smoke_split_name=str(values["smoke_split_name"]),
        smoke_fold_id=str(values["smoke_fold_id"]),
        hidden_dim=int(values["hidden_dim"]),
        capacity_matched_hidden_dim=int(values["capacity_matched_hidden_dim"]),
        learning_rate=float(values["learning_rate"]),
        weight_decay=float(values["weight_decay"]),
        epochs=int(values["epochs"]),
        batch_size=int(values["batch_size"]),
        tie_tolerance=float(values["tie_tolerance"]),
    )


def sample_protocol_pairs(
    benchmark: pd.DataFrame,
    train_proteins: list[str],
    protocol: str,
    config: V05Config,
) -> pd.DataFrame:
    if protocol == "S512":
        return sample_rank_pairs(
            benchmark,
            train_proteins,
            pairs_per_protein=512,
            seed=config.seed,
            tie_tolerance=config.tie_tolerance,
        )
    pair_count = {"D4096": 4096, "D16384": 16384}[protocol]
    return sample_dense_pairs(
        benchmark,
        train_proteins,
        pairs_per_protein=pair_count,
        seed=config.seed,
        protocol=protocol,
        tie_tolerance=config.tie_tolerance,
    )


def _model_factory(label: str, protein_dim: int, config: V05Config):
    if label == "M0":
        return CandidateDNAOnly(config.hidden_dim), CandidateDNAOnly.model_name
    if label == "M1c":
        return (
            ProteinCandidate(protein_dim, config.capacity_matched_hidden_dim),
            ProteinCandidate.model_name,
        )
    if label == "M2":
        return TargetCandidateOnly(config.hidden_dim), TargetCandidateOnly.model_name
    if label == "M3":
        return ProteinTargetCandidate(protein_dim, config.hidden_dim), ProteinTargetCandidate.model_name
    raise ValueError(f"Unknown dense model label {label}")


def _metric_for_prediction(
    benchmark: pd.DataFrame,
    protein: str,
    prediction: np.ndarray,
) -> dict[str, float | int]:
    group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
    metrics = compute_ranking_metrics(
        pd.DataFrame(
            {
                "experimental_score": group["experimental_score"],
                "prediction": prediction,
                "canonical_7mer": group["candidate_dna"],
            }
        ),
        "experimental_score",
        "prediction",
    )
    return {
        "spearman": metrics.spearman,
        "ndcg_1pct": metrics.ndcg_1pct,
        "ndcg_5pct": metrics.ndcg_5pct,
        "pairwise_accuracy": metrics.pairwise_accuracy,
        "top1pct_recovery": metrics.top1pct_recovery,
        "n_rc_units": metrics.n_rc_classes,
    }


def _corr(first: np.ndarray, second: np.ndarray) -> float:
    if np.std(first) == 0.0 or np.std(second) == 0.0:
        return np.nan
    return float(np.corrcoef(first, second)[0, 1])


def _shuffled_predictions(
    model,
    model_name: str,
    label: str,
    protocol: str,
    cache: dict[str, dict[str, Any]],
    test_protein: str,
    shuffled_protein: str,
) -> list[dict[str, Any]]:
    entry = cache[test_protein]
    original = predict_model(model, model_name, entry)
    rows: list[dict[str, Any]] = []
    if label == "M0":
        return [
            {
                "model": "M0",
                "protocol": protocol,
                "dbp_id": test_protein,
                "shuffle_type": "none",
                "permutation_source": "",
                "prediction_correlation": np.nan,
                "mean_abs_score_change": 0.0,
                "retrained": False,
            }
        ]
    if label in {"M1c", "M3"}:
        protein_shuffled = dict(entry)
        protein_shuffled["embedding"] = cache[shuffled_protein]["embedding"]
        shuffled = predict_model(model, model_name, protein_shuffled)
        rows.append(
            {
                "model": label,
                "protocol": protocol,
                "dbp_id": test_protein,
                "shuffle_type": "protein",
                "permutation_source": shuffled_protein,
                "prediction_correlation": _corr(original, shuffled),
                "mean_abs_score_change": float(np.mean(np.abs(original - shuffled))),
                "retrained": False,
            }
        )
    if label in {"M2", "M3"}:
        target_shuffled = dict(entry)
        target_shuffled["target"] = cache[shuffled_protein]["target"]
        target_shuffled["target_features"] = cache[shuffled_protein]["target_features"]
        shuffled_target = predict_model(model, model_name, target_shuffled)
        rows.append(
            {
                "model": label,
                "protocol": protocol,
                "dbp_id": test_protein,
                "shuffle_type": "target",
                "permutation_source": shuffled_protein,
                "prediction_correlation": _corr(original, shuffled_target),
                "mean_abs_score_change": float(np.mean(np.abs(original - shuffled_target))),
                "retrained": False,
            }
        )
    return rows


def run_dense_smoke(config: V05Config | None = None) -> dict[str, pd.DataFrame]:
    config = dense_config_from_registered() if config is None else config
    set_seed(config.seed)
    benchmark, targets, embedding_map, splits = load_v05_data()
    fold_id, train_proteins, test_proteins = select_untouched_smoke_fold(
        splits,
        still_untouched={"DBP48", "DBP6", "DBP9"},
        split_name=config.smoke_split_name,
    )
    if fold_id != config.smoke_fold_id:
        raise AssertionError(f"Registered fold {config.smoke_fold_id} disagrees with manifest selection {fold_id}")
    cache = build_cache(benchmark, targets, embedding_map, sorted(set(train_proteins + test_proteins)))
    protein_dim = len(next(iter(embedding_map.values())))
    all_eval_rows: list[dict[str, Any]] = []
    health_rows: list[dict[str, Any]] = []
    shuffle_rows: list[dict[str, Any]] = []
    prediction_store: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    history_rows: list[pd.DataFrame] = []
    runtime_rows: list[dict[str, Any]] = []
    pair_rows: list[pd.DataFrame] = []

    for protocol in ("S512", "D4096", "D16384"):
        pairs = sample_protocol_pairs(benchmark, train_proteins, protocol, config)
        pair_rows.append(pairs.assign(protocol=protocol))
        protocol_config = replace(
            config,
            pair_count_per_protein={"S512": 512, "D4096": 4096, "D16384": 16384}[protocol],
        )
        for label in ("M0", "M1c", "M2", "M3"):
            set_seed(config.seed)
            model, canonical_name = _model_factory(label, protein_dim, config)
            model, history, runtime = train_model(
                model,
                canonical_name,
                cache,
                pairs,
                train_proteins,
                protocol_config,
            )
            history = history.copy()
            history["reported_model"] = label
            history["protocol"] = protocol
            history_rows.append(history)
            runtime_rows.append(
                {
                    "model": label,
                    "protocol": protocol,
                    "runtime_seconds": runtime,
                    "training_proteins": "|".join(train_proteins),
                    "test_proteins": "|".join(test_proteins),
                }
            )
            train_predictions = {
                protein: predict_model(model, canonical_name, cache[protein]) for protein in train_proteins
            }
            test_predictions = {
                protein: predict_model(model, canonical_name, cache[protein]) for protein in test_proteins
            }
            prediction_store[(protocol, label)] = test_predictions
            for protein in train_proteins:
                metrics = _metric_for_prediction(benchmark, protein, train_predictions[protein])
                all_eval_rows.append(
                    {
                        "model": label,
                        "protocol": protocol,
                        "dbp_id": protein,
                        "partition": "train",
                        **metrics,
                    }
                )
            for protein in test_proteins:
                metrics = _metric_for_prediction(benchmark, protein, test_predictions[protein])
                all_eval_rows.append(
                    {
                        "model": label,
                        "protocol": protocol,
                        "dbp_id": protein,
                        "partition": "test",
                        **metrics,
                    }
                )
            train_rhos = [
                row["spearman"]
                for row in all_eval_rows
                if row["model"] == label and row["protocol"] == protocol and row["partition"] == "train"
            ]
            test_values = np.concatenate(list(test_predictions.values()))
            last_loss = float(history["mean_pairwise_loss"].iloc[-1])
            first_loss = float(history["mean_pairwise_loss"].iloc[0])
            steps_per_epoch = sum(
                int(np.ceil(len(pairs.loc[pairs["protein_id"].eq(protein)]) / config.batch_size))
                for protein in train_proteins
            )
            health_rows.append(
                {
                    "model": label,
                    "protocol": protocol,
                    "training_pair_count": len(pairs),
                    "optimizer_steps": steps_per_epoch * config.epochs,
                    "first_epoch_loss": first_loss,
                    "final_epoch_loss": last_loss,
                    "loss_reduction": first_loss - last_loss,
                    "train_macro_median_spearman": float(np.nanmedian(train_rhos)),
                    "test_macro_median_spearman": float(
                        np.nanmedian(
                            [
                                row["spearman"]
                                for row in all_eval_rows
                                if row["model"] == label
                                and row["protocol"] == protocol
                                and row["partition"] == "test"
                            ]
                        )
                    ),
                    "test_prediction_variance": float(np.var(test_values)),
                    "nan_inf_count": int(np.sum(~np.isfinite(test_values))),
                    "runtime_seconds": runtime,
                    "status": "complete",
                }
            )
            shuffled_source = sorted(set(train_proteins + test_proteins))
            source = shuffled_source[(shuffled_source.index(test_proteins[0]) + 1) % len(shuffled_source)]
            shuffle_rows.extend(
                _shuffled_predictions(model, canonical_name, label, protocol, cache, test_proteins[0], source)
            )

    evaluation = pd.DataFrame(all_eval_rows)
    health = pd.DataFrame(health_rows)
    history = pd.concat(history_rows, ignore_index=True)
    pairs = pd.concat(pair_rows, ignore_index=True)
    shuffle = pd.DataFrame(shuffle_rows)
    delta_rows = []
    for label in ("M0", "M1c", "M2", "M3"):
        for protein in test_proteins:
            values = evaluation.loc[
                evaluation["model"].eq(label)
                & evaluation["dbp_id"].eq(protein)
                & evaluation["partition"].eq("test"),
                ["protocol", "spearman"],
            ].set_index("protocol")["spearman"]
            delta_rows.append(
                {
                    "model": label,
                    "dbp_id": protein,
                    "S512_test_rho": float(values["S512"]),
                    "D4096_test_rho": float(values["D4096"]),
                    "D16384_test_rho": float(values["D16384"]),
                    "D4096_minus_S512": float(values["D4096"] - values["S512"]),
                    "D16384_minus_S512": float(values["D16384"] - values["S512"]),
                }
            )
    delta = pd.DataFrame(delta_rows)
    coverage = pd.read_csv(OUTPUT_DIR / "pair_sampling_audit.csv")
    chain_rows = []
    for protocol in ("S512", "D4096", "D16384"):
        protocol_cov = float(
            coverage.loc[coverage["protocol"].eq(protocol), "candidate_coverage"].median()
        )
        for label in ("M1c", "M2", "M3"):
            train_rho = float(
                health.loc[health["model"].eq(label) & health["protocol"].eq(protocol), "train_macro_median_spearman"].iloc[0]
            )
            test_rho = float(
                health.loc[health["model"].eq(label) & health["protocol"].eq(protocol), "test_macro_median_spearman"].iloc[0]
            )
            sensitivity = shuffle.loc[
                shuffle["model"].eq(label) & shuffle["protocol"].eq(protocol), "prediction_correlation"
            ].dropna()
            chain_rows.append(
                {
                    "model": label,
                    "protocol": protocol,
                    "median_candidate_coverage": protocol_cov,
                    "train_rho": train_rho,
                    "test_rho": test_rho,
                    "median_shuffle_prediction_correlation": float(sensitivity.median())
                    if not sensitivity.empty
                    else np.nan,
                }
            )
    chain = pd.DataFrame(chain_rows)
    info = pd.DataFrame(
        [
            {
                "split_name": config.smoke_split_name,
                "fold_id": fold_id,
                "train_proteins": "|".join(train_proteins),
                "test_proteins": "|".join(test_proteins),
                "seed": config.seed,
                "status": "NOT PRIMARY SCIENTIFIC RESULT",
                "exposure_status": "development-exposed-for-dense-supervision",
            }
        ]
    )
    return {
        "evaluation": evaluation,
        "health": health,
        "history": history,
        "runtime": pd.DataFrame(runtime_rows),
        "pairs": pairs,
        "shuffle": shuffle,
        "delta": delta,
        "chain": chain,
        "info": info,
    }


def write_dense_outputs(outputs: dict[str, pd.DataFrame]) -> None:
    ensure_dir(OUTPUT_DIR)
    ensure_dir(LOG_DIR)
    outputs["evaluation"].to_csv(OUTPUT_DIR / "smoke_results.csv", index=False)
    outputs["health"].to_csv(OUTPUT_DIR / "training_health.csv", index=False)
    outputs["delta"].to_csv(OUTPUT_DIR / "supervision_delta_summary.csv", index=False)
    outputs["shuffle"].to_csv(OUTPUT_DIR / "shuffle_diagnostics.csv", index=False)
    outputs["chain"].to_csv(OUTPUT_DIR / "mechanism_chain.csv", index=False)
    outputs["info"].to_csv(OUTPUT_DIR / "smoke_split_info.csv", index=False)
    outputs["history"].to_csv(OUTPUT_DIR / "smoke_loss_history.csv", index=False)
    outputs["runtime"].to_csv(OUTPUT_DIR / "smoke_runtime.csv", index=False)
    outputs["pairs"].to_csv(OUTPUT_DIR / "smoke_training_pairs.csv", index=False)
    log_parts = [
        "Phase 6A dense supervision smoke",
        outputs["info"].to_string(index=False),
        outputs["health"].to_string(index=False),
        outputs["delta"].to_string(index=False),
        outputs["shuffle"].to_string(index=False),
        outputs["chain"].to_string(index=False),
    ]
    (LOG_DIR / "dense_supervision_smoke.log").write_text("\n\n".join(log_parts), encoding="utf-8")
