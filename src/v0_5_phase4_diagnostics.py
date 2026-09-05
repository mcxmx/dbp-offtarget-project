from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

from src.sequence_equivalence import canonical_rc
from src.utils import project_root, sequence_identity
from src.v0_4_evaluation import compute_ranking_metrics
from src.v0_5_models import (
    CandidateDNAOnly,
    ProteinCandidate,
    ProteinTargetCandidate,
    TargetCandidateOnly,
    score_rc,
    target_edit_control,
    target_hamming_control,
    target_kmer_overlap_control,
)
from src.v0_5_primary_evaluation import (
    MODEL_ORDER,
    canonical_model_name,
    get_fold_partitions,
    model_factory,
)
from src.v0_5_training import (
    V05Config,
    build_cache,
    load_v05_data,
    predict_model,
    sample_rank_pairs,
    set_seed,
    train_model,
)


ROOT = project_root()
DESIGNED_IDS = ("DBP1", "DBP3", "DBP35", "DBP48", "DBP5", "DBP6", "DBP9")
PRIMARY_RESULT_FILES = (
    "results/v0_5/primary_seed_level_results.csv",
    "results/v0_5/primary_per_protein_results.csv",
    "results/v0_5/primary_macro_summary.csv",
    "results/v0_5/strict_component_seed_level_results.csv",
    "results/v0_5/strict_component_per_protein_results.csv",
    "results/v0_5/strict_component_macro_summary.csv",
    "results/v0_5/target_relative_controls_full.csv",
    "results/v0_5/baseline_context_table.csv",
    "results/v0_5/primary_training_health.csv",
    "results/v0_5/V0_5_PRIMARY_RESULTS.md",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_primary_frozen_manifest(
    root: Path = ROOT,
    *,
    source_commit: str,
    primary_result_commit: str,
    tag: str,
    seeds: Iterable[int],
    config_path: str = "metadata/v0_5/v0_5_model_config.json",
    split_manifest: str = "metadata/v0_5/v0_5_split_manifest.csv",
) -> Path:
    manifest = root / "results" / "v0_5" / "PRIMARY_RESULTS_FROZEN_MANIFEST.txt"
    lines = [
        "# v0.5 primary result freeze",
        f"source_commit={source_commit}",
        f"primary_result_commit={primary_result_commit}",
        f"freeze_tag={tag}",
        f"seeds={'|'.join(str(seed) for seed in seeds)}",
        f"config={config_path}",
        f"split_manifest={split_manifest}",
        "hash_algorithm=sha256",
        "",
        "# path sha256",
    ]
    for relative in PRIMARY_RESULT_FILES:
        path = root / relative
        if not path.exists():
            raise FileNotFoundError(path)
        lines.append(f"{relative} {sha256_file(path)}")
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def safe_spearman(left: Iterable[float], right: Iterable[float]) -> float:
    left_array = np.asarray(list(left), dtype=float)
    right_array = np.asarray(list(right), dtype=float)
    valid = np.isfinite(left_array) & np.isfinite(right_array)
    left_array = left_array[valid]
    right_array = right_array[valid]
    if len(left_array) < 2 or np.unique(left_array).size < 2 or np.unique(right_array).size < 2:
        return np.nan
    return float(spearmanr(left_array, right_array).statistic)


def safe_pearson(left: Iterable[float], right: Iterable[float]) -> float:
    left_array = np.asarray(list(left), dtype=float)
    right_array = np.asarray(list(right), dtype=float)
    valid = np.isfinite(left_array) & np.isfinite(right_array)
    left_array = left_array[valid]
    right_array = right_array[valid]
    if len(left_array) < 2 or np.std(left_array) == 0 or np.std(right_array) == 0:
        return np.nan
    return float(np.corrcoef(left_array, right_array)[0, 1])


def deterministic_protein_permutation(protein_ids: Iterable[str]) -> dict[str, str]:
    ordered = sorted(set(protein_ids))
    if len(ordered) < 2:
        raise ValueError("At least two proteins are required for a permutation")
    return {
        protein: ordered[(index + 1) % len(ordered)]
        for index, protein in enumerate(ordered)
    }


def deterministic_target_permutation(protein_ids: Iterable[str]) -> dict[str, str]:
    ordered = [protein for protein in DESIGNED_IDS if protein in set(protein_ids)]
    if len(ordered) < 2:
        raise ValueError("At least two designed proteins are required for a target permutation")
    return {
        protein: ordered[(index + 1) % len(ordered)]
        for index, protein in enumerate(ordered)
    }


def _prediction_with_overrides(
    model: torch.nn.Module,
    model_name: str,
    entry: dict[str, Any],
    *,
    protein_embedding: np.ndarray,
    target: str,
    chunk_size: int = 512,
) -> np.ndarray:
    output: list[np.ndarray] = []
    candidates = entry["candidate_dna"]
    embedding = torch.from_numpy(np.asarray(protein_embedding, dtype=np.float32))
    for start in range(0, len(candidates), chunk_size):
        output.append(
            score_rc(
                model,
                model_name,
                candidates[start : start + chunk_size],
                protein_embedding=embedding,
                target=target,
            )
        )
    return np.concatenate(output)


def _pair_sampling_coverage(
    benchmark: pd.DataFrame,
    pairs: pd.DataFrame,
    *,
    split_type: str,
    fold_id: str,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for protein, group_pairs in pairs.groupby("protein_id", sort=True):
        group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
        indices = np.unique(
            np.concatenate(
                [
                    group_pairs["left_index"].to_numpy(dtype=int),
                    group_pairs["right_index"].to_numpy(dtype=int),
                ]
            )
        )
        ranks = pd.Series(group["experimental_score"]).rank(method="average", pct=True).to_numpy()
        used_ranks = ranks[indices]
        rank_bins = np.floor(np.minimum(used_ranks, 1.0 - 1e-12) * 10).astype(int)
        rows.append(
            {
                "split_type": split_type,
                "fold_id": fold_id,
                "seed": int(seed),
                "protein_id": protein,
                "pairs_sampled": int(len(group_pairs)),
                "n_rc_classes": int(len(group)),
                "unique_candidates_in_pairs": int(len(indices)),
                "candidate_coverage_fraction": float(len(indices) / len(group)),
                "min_rank_percentile": float(np.min(used_ranks)),
                "max_rank_percentile": float(np.max(used_ranks)),
                "q05_rank_percentile": float(np.quantile(used_ranks, 0.05)),
                "q95_rank_percentile": float(np.quantile(used_ranks, 0.95)),
                "n_decile_bins_covered": int(np.unique(rank_bins).size),
                "label_quantile_min": float(np.min(used_ranks)),
                "label_quantile_max": float(np.max(used_ranks)),
                "sampling_protocol": "frozen 40% easy, 35% medium, 25% hard within-protein pairs",
            }
        )
    return pd.DataFrame(rows)


def replay_primary_predictions(
    config: V05Config,
    *,
    split_name: str = "protein_cluster_loco",
    seeds: Iterable[int] | None = None,
) -> dict[str, pd.DataFrame]:
    """Replay the frozen primary training contract for candidate-level diagnostics.

    This intentionally replays the fixed primary runs without changing any
    architecture, hyperparameter, split, seed, or sampling rule. It is not a
    new model evaluation and its outputs are kept outside the frozen primary
    result files.
    """
    benchmark, targets, embedding_map, splits = load_v05_data()
    seeds = tuple(config.evaluation_seeds if seeds is None else seeds)
    protein_dim = len(next(iter(embedding_map.values())))
    prediction_frames: list[pd.DataFrame] = []
    train_test_rows: list[dict[str, Any]] = []
    health_rows: list[dict[str, Any]] = []
    pair_coverage_rows: list[pd.DataFrame] = []
    shuffled_protein_rows: list[dict[str, Any]] = []
    shuffled_target_rows: list[dict[str, Any]] = []
    target_map = targets.set_index("dbp_id")["primary_target"].to_dict()
    protein_permutation = deterministic_protein_permutation(DESIGNED_IDS)
    target_permutation = deterministic_target_permutation(DESIGNED_IDS)

    fold_ids = sorted(splits.loc[splits["split_name"].eq(split_name), "fold_id"].unique())
    for fold_id in fold_ids:
        train_proteins, test_proteins = get_fold_partitions(splits, split_name, fold_id)
        all_fold_proteins = sorted(set(train_proteins + test_proteins))
        cache = build_cache(benchmark, targets, embedding_map, all_fold_proteins)
        fold_prediction_frames: list[pd.DataFrame] = []
        for seed in seeds:
            run_config = replace(config, seed=int(seed))
            pairs = sample_rank_pairs(
                benchmark,
                train_proteins,
                pairs_per_protein=config.pair_count_per_protein,
                seed=int(seed),
                tie_tolerance=config.tie_tolerance,
            )
            pair_coverage_rows.append(
                _pair_sampling_coverage(
                    benchmark,
                    pairs,
                    split_type=split_name,
                    fold_id=fold_id,
                    seed=int(seed),
                )
            )
            for short_name in MODEL_ORDER:
                set_seed(int(seed))
                model = model_factory(short_name, protein_dim, run_config)
                model, history, runtime = train_model(
                    model,
                    canonical_model_name(short_name),
                    cache,
                    pairs,
                    train_proteins,
                    run_config,
                )
                canonical_name = canonical_model_name(short_name)
                test_scores: dict[str, np.ndarray] = {}
                for partition, proteins in (("train", train_proteins), ("test", test_proteins)):
                    for protein in proteins:
                        scores = predict_model(model, canonical_name, cache[protein])
                        truth = cache[protein]["truth"]
                        metrics = compute_ranking_metrics(
                            pd.DataFrame({"truth": truth, "prediction": scores}),
                            "truth",
                            "prediction",
                        )
                        train_test_rows.append(
                            {
                                "split_type": split_name,
                                "fold_id": fold_id,
                                "seed": int(seed),
                                "partition": partition,
                                "dbp_id": protein,
                                "model": short_name,
                                "spearman": metrics.spearman,
                                "n_rc_classes": metrics.n_rc_classes,
                                "diagnostic_only": True,
                                "used_for_model_selection": False,
                            }
                        )
                        if partition == "test":
                            test_scores[protein] = scores
                            fold_prediction_frames.append(
                                pd.DataFrame(
                                    {
                                        "split_type": split_name,
                                        "fold_id": fold_id,
                                        "seed": int(seed),
                                        "partition": "test",
                                        "dbp_id": protein,
                                        "canonical_7mer": cache[protein]["candidate_dna"],
                                        "prediction_model": short_name,
                                        "prediction_score": scores,
                                        "replay_status": "exact_frozen_primary_protocol",
                                    }
                                )
                            )
                all_predictions = np.concatenate(list(test_scores.values()))
                health_rows.append(
                    {
                        "split_type": split_name,
                        "fold_id": fold_id,
                        "seed": int(seed),
                        "model": short_name,
                        "training_proteins": "|".join(train_proteins),
                        "test_proteins": "|".join(test_proteins),
                        "training_pair_count": int(len(pairs)),
                        "runtime_seconds": float(runtime),
                        "first_epoch_loss": float(history["mean_pairwise_loss"].iloc[0]),
                        "last_epoch_loss": float(history["mean_pairwise_loss"].iloc[-1]),
                        "prediction_variance": float(np.var(all_predictions)),
                        "nan_inf_count": int(np.sum(~np.isfinite(all_predictions))),
                        "diagnostic_replay": True,
                    }
                )
                if short_name == "M3":
                    for protein in test_proteins:
                        entry = cache[protein]
                        original = test_scores[protein]
                        shuffled_source = protein_permutation[protein]
                        shuffled_p = _prediction_with_overrides(
                            model,
                            canonical_name,
                            entry,
                            protein_embedding=embedding_map[shuffled_source],
                            target=entry["target"],
                        )
                        shuffled_target_source = target_permutation[protein]
                        shuffled_t = _prediction_with_overrides(
                            model,
                            canonical_name,
                            entry,
                            protein_embedding=embedding_map[protein],
                            target=target_map[shuffled_target_source],
                        )
                        for condition, scores, source in (
                            ("original", original, protein),
                            ("shuffled_protein", shuffled_p, shuffled_source),
                        ):
                            shuffled_protein_rows.append(
                                {
                                    "split_type": split_name,
                                    "fold_id": fold_id,
                                    "seed": int(seed),
                                    "dbp_id": protein,
                                    "condition": condition,
                                    "replacement_protein_id": source,
                                    "prediction_correlation_to_original": (
                                        1.0 if condition == "original" else safe_pearson(original, scores)
                                    ),
                                    "experimental_spearman": safe_spearman(entry["truth"], scores),
                                    "mean_absolute_score_change": float(np.mean(np.abs(scores - original))),
                                    "prediction_variance": float(np.var(scores)),
                                    "retrained": False,
                                    "used_for_model_selection": False,
                                }
                            )
                        for condition, scores, source in (
                            ("original", original, protein),
                            ("shuffled_target", shuffled_t, shuffled_target_source),
                        ):
                            shuffled_target_rows.append(
                                {
                                    "split_type": split_name,
                                    "fold_id": fold_id,
                                    "seed": int(seed),
                                    "dbp_id": protein,
                                    "condition": condition,
                                    "replacement_target_dbp_id": source,
                                    "prediction_correlation_to_original": (
                                        1.0 if condition == "original" else safe_pearson(original, scores)
                                    ),
                                    "experimental_spearman": safe_spearman(entry["truth"], scores),
                                    "mean_absolute_score_change": float(np.mean(np.abs(scores - original))),
                                    "prediction_variance": float(np.var(scores)),
                                    "retrained": False,
                                    "used_for_model_selection": False,
                                }
                            )
        prediction_frames.extend(fold_prediction_frames)

    return {
        "predictions": pd.concat(prediction_frames, ignore_index=True),
        "train_test": pd.DataFrame(train_test_rows),
        "health": pd.DataFrame(health_rows),
        "pair_coverage": pd.concat(pair_coverage_rows, ignore_index=True),
        "shuffled_protein": pd.DataFrame(shuffled_protein_rows),
        "shuffled_target": pd.DataFrame(shuffled_target_rows),
    }


def build_wide_predictions(
    replay_predictions: pd.DataFrame,
    benchmark: pd.DataFrame,
) -> pd.DataFrame:
    wide = (
        replay_predictions.groupby(
            ["dbp_id", "canonical_7mer", "prediction_model"], as_index=False
        )["prediction_score"]
        .mean()
        .pivot_table(
            index=["dbp_id", "canonical_7mer"],
            columns="prediction_model",
            values="prediction_score",
            aggfunc="first",
        )
        .reset_index()
    )
    wide.columns.name = None
    wide = wide.rename(columns={model: f"{model}_score" for model in MODEL_ORDER})
    truth = benchmark.rename(
        columns={"protein_id": "dbp_id", "candidate_dna": "canonical_7mer"}
    )[["dbp_id", "canonical_7mer", "experimental_score"]]
    return truth.merge(wide, on=["dbp_id", "canonical_7mer"], validate="one_to_one")


def add_target_controls(
    wide: pd.DataFrame,
    targets: pd.DataFrame,
) -> pd.DataFrame:
    target_map = targets.set_index("dbp_id")["primary_target"].to_dict()
    output = wide.copy()
    values: dict[str, list[float]] = {
        "target_hamming": [],
        "target_edit": [],
        "target_kmer_overlap": [],
    }
    for _, row in output.iterrows():
        target = target_map[row["dbp_id"]]
        candidate = row["canonical_7mer"]
        values["target_hamming"].append(target_hamming_control(target, candidate))
        values["target_edit"].append(target_edit_control(target, candidate))
        values["target_kmer_overlap"].append(target_kmer_overlap_control(target, candidate))
    for name, series in values.items():
        output[name] = series
    output["experimental_percentile"] = output.groupby("dbp_id")["experimental_score"].rank(pct=True)
    for model in MODEL_ORDER:
        score_col = f"{model}_score"
        output[f"{model}_percentile"] = output.groupby("dbp_id")[score_col].rank(pct=True)
    return output


def load_reference_hard_cases(root: Path = ROOT) -> pd.DataFrame:
    path = root / "results" / "v0_4_2" / "tables" / "baseline_failure_cases_deeppbs_completed_v0_4_2.parquet"
    reference = pd.read_parquet(path)[
        ["protein_id", "canonical_7mer", "failure_category"]
    ].rename(columns={"failure_category": "existing_hard_case_category"})
    if len(reference) != 1515 or reference.duplicated(["protein_id", "canonical_7mer"]).any():
        raise AssertionError("The v0.3.1 disagreement reference set is not the expected 1,515 unique RC-class rows")
    return reference


def build_hard_case_tables(
    wide: pd.DataFrame,
    reference: pd.DataFrame,
    output_dir: Path,
) -> dict[str, pd.DataFrame]:
    reference = reference.rename(columns={"protein_id": "dbp_id"})
    hard = reference.merge(wide, on=["dbp_id", "canonical_7mer"], how="left")
    hard = hard.rename(columns={"canonical_7mer": "candidate"})
    for model in MODEL_ORDER:
        hard[f"{model}_resolved"] = hard[f"{model}_percentile"] >= 0.90
        hard[f"{model}_evaluable"] = hard[f"{model}_score"].notna()
    hard["rc_class"] = hard["candidate"].map(canonical_rc)
    hard["experimental_E_score"] = hard["experimental_score"]
    hard["target_similarity_kmer_overlap"] = hard["target_kmer_overlap"]
    base_columns = [
        "dbp_id",
        "candidate",
        "rc_class",
        "experimental_E_score",
        "M0_score",
        "M1_score",
        "M1c_score",
        "M2_score",
        "M3_score",
        "target_hamming",
        "target_edit",
        "target_kmer_overlap",
        "existing_hard_case_category",
    ]
    subset_masks = {
        "m1_fail_m3_success": (~hard["M1_resolved"]) & hard["M3_resolved"],
        "m1c_fail_m3_success": (~hard["M1c_resolved"]) & hard["M3_resolved"],
        "m2_fail_m3_success": (~hard["M2_resolved"]) & hard["M3_resolved"],
        "joint_controls_fail_m3_success": (
            (~hard["M1c_resolved"]) & (~hard["M2_resolved"]) & hard["M3_resolved"]
        ),
        "all_current_models_fail": ~hard[[f"{model}_resolved" for model in MODEL_ORDER]].any(axis=1),
    }
    subset_tables: dict[str, pd.DataFrame] = {}
    for name, mask in subset_masks.items():
        table = hard.loc[mask, base_columns].copy().sort_values(["dbp_id", "candidate"])
        table.to_csv(output_dir / f"{name}.csv", index=False)
        subset_tables[name] = table
    summary_rows: list[dict[str, Any]] = []
    for name, table in subset_tables.items():
        summary_rows.append(
            {
                "subset": name,
                "count": int(len(table)),
                "fraction_of_reference": float(len(table) / len(hard)),
                "proteins_represented": "|".join(sorted(table["dbp_id"].unique())),
                "n_proteins": int(table["dbp_id"].nunique()),
                "experimental_E_score_median": float(table["experimental_E_score"].median()) if len(table) else np.nan,
                "target_hamming_median": float(table["target_hamming"].median()) if len(table) else np.nan,
                "target_edit_median": float(table["target_edit"].median()) if len(table) else np.nan,
                "target_kmer_overlap_median": float(table["target_kmer_overlap"].median()) if len(table) else np.nan,
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output_dir.parent / "hard_case_subset_summary.csv", index=False)
    subset_by_protein = (
        hard.assign(
            m1_fail_m3_success=subset_masks["m1_fail_m3_success"],
            m1c_fail_m3_success=subset_masks["m1c_fail_m3_success"],
            m2_fail_m3_success=subset_masks["m2_fail_m3_success"],
            joint_controls_fail_m3_success=subset_masks["joint_controls_fail_m3_success"],
            all_current_models_fail=subset_masks["all_current_models_fail"],
        )
        .groupby("dbp_id", as_index=False)[
            [
                "m1_fail_m3_success",
                "m1c_fail_m3_success",
                "m2_fail_m3_success",
                "joint_controls_fail_m3_success",
                "all_current_models_fail",
            ]
        ]
        .sum()
    )
    subset_by_protein.to_csv(output_dir.parent / "hard_case_subset_by_protein.csv", index=False)
    by_protein_rows: list[dict[str, Any]] = []
    for protein, group in hard.groupby("dbp_id", sort=True):
        row: dict[str, Any] = {"dbp_id": protein, "hard_cases": int(len(group))}
        for model in MODEL_ORDER:
            row[model] = int(group[f"{model}_resolved"].sum())
        by_protein_rows.append(row)
    by_protein = pd.DataFrame(by_protein_rows)
    by_protein.to_csv(output_dir.parent / "hard_case_resolution_by_protein.csv", index=False)
    return {
        "hard": hard,
        "summary": summary,
        "by_protein": by_protein,
        **subset_tables,
    }


def build_hard_case_model_summary(hard: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    rows = []
    for model in MODEL_ORDER:
        evaluable = hard[f"{model}_evaluable"]
        resolved = hard[f"{model}_resolved"] & evaluable
        rows.append(
            {
                "model": model,
                "eligible": int(evaluable.sum()),
                "resolved": int(resolved.sum()),
                "unresolved": int((evaluable & ~resolved).sum()),
                "not_evaluable": int((~evaluable).sum()),
                "resolution_rate": float(resolved.sum() / evaluable.sum()) if evaluable.sum() else np.nan,
                "reference_total": int(len(hard)),
                "resolution_definition": "existing v0.3.1 hard cases; prediction percentile >= 0.90 after seed-mean aggregation",
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)
    return result


def target_similarity_correlations(wide: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for protein, group in wide.groupby("dbp_id", sort=True):
        for model in MODEL_ORDER:
            rows.append(
                {
                    "dbp_id": protein,
                    "model": model,
                    "spearman_prediction_vs_target_kmer_overlap": safe_spearman(
                        group[f"{model}_score"], group["target_kmer_overlap"]
                    ),
                    "spearman_prediction_vs_experimental": safe_spearman(
                        group[f"{model}_score"], group["experimental_score"]
                    ),
                    "n_rc_classes": int(len(group)),
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)
    return result


def target_similarity_bin_performance(wide: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for protein, group in wide.groupby("dbp_id", sort=True):
        q33, q67 = np.quantile(group["target_kmer_overlap"], [1 / 3, 2 / 3])
        work = group.copy()
        work["target_similarity_bin"] = np.select(
            [
                work["target_kmer_overlap"] <= q33,
                work["target_kmer_overlap"] <= q67,
            ],
            ["low", "medium"],
            default="high",
        )
        for bin_name, bin_group in work.groupby("target_similarity_bin", sort=True):
            for model in ("M1c", "M2", "M3"):
                rows.append(
                    {
                        "dbp_id": protein,
                        "target_similarity_bin": bin_name,
                        "bin_lower_threshold": float(q33),
                        "bin_upper_threshold": float(q67),
                        "model": model,
                        "spearman": safe_spearman(
                            bin_group[f"{model}_score"], bin_group["experimental_score"]
                        ),
                        "n_rc_classes": int(len(bin_group)),
                        "bin_definition": "per-protein target-kmer-overlap tertiles fixed before outcome inspection",
                    }
                )
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)
    return result


def partial_spearman_target_diagnostic(wide: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for protein, group in wide.groupby("dbp_id", sort=True):
        control_rank = pd.Series(group["target_kmer_overlap"]).rank(method="average").to_numpy()
        experimental_rank = pd.Series(group["experimental_score"]).rank(method="average").to_numpy()
        target_slope, target_intercept = np.polyfit(control_rank, experimental_rank, 1)
        experimental_residual = experimental_rank - (target_slope * control_rank + target_intercept)
        for model in MODEL_ORDER:
            prediction_rank = pd.Series(group[f"{model}_score"]).rank(method="average").to_numpy()
            pred_slope, pred_intercept = np.polyfit(control_rank, prediction_rank, 1)
            prediction_residual = prediction_rank - (pred_slope * control_rank + pred_intercept)
            rows.append(
                {
                    "dbp_id": protein,
                    "model": model,
                    "partial_spearman_residual_association": safe_spearman(
                        experimental_residual, prediction_residual
                    ),
                    "ordinary_spearman": safe_spearman(
                        group["experimental_score"], group[f"{model}_score"]
                    ),
                    "control": "TargetKmerOverlap",
                    "diagnostic_status": "SECONDARY_DIAGNOSTIC",
                    "n_rc_classes": int(len(group)),
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)
    return result


def embedding_similarity_table(
    embedding_map: dict[str, np.ndarray],
    targets: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    sequence_map = targets.set_index("dbp_id")["protein_sequence"].to_dict()
    rows = []
    for left in DESIGNED_IDS:
        for right in DESIGNED_IDS:
            left_vector = embedding_map[left]
            right_vector = embedding_map[right]
            denominator = np.linalg.norm(left_vector) * np.linalg.norm(right_vector)
            rows.append(
                {
                    "dbp_id_i": left,
                    "dbp_id_j": right,
                    "embedding_cosine_similarity": float(np.dot(left_vector, right_vector) / denominator),
                    "embedding_euclidean_distance": float(np.linalg.norm(left_vector - right_vector)),
                    "sequence_identity": float(sequence_identity(sequence_map[left], sequence_map[right])),
                    "same_dbp": left == right,
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)
    return result


def embedding_pca_table(
    embedding_map: dict[str, np.ndarray],
    targets: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    ordered = [protein for protein in DESIGNED_IDS if protein in embedding_map]
    matrix = np.vstack([embedding_map[protein] for protein in ordered]).astype(float)
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    coordinates = centered @ vt[:2].T
    variance = singular_values**2
    variance = variance / variance.sum() if variance.sum() else variance
    cluster_map = targets.set_index("dbp_id")["protein_cluster"].to_dict()
    result = pd.DataFrame(
        {
            "dbp_id": ordered,
            "protein_cluster": [cluster_map[protein] for protein in ordered],
            "pc1": coordinates[:, 0],
            "pc2": coordinates[:, 1],
            "pc1_explained_fraction": float(variance[0]) if len(variance) > 0 else np.nan,
            "pc2_explained_fraction": float(variance[1]) if len(variance) > 1 else np.nan,
            "embedding_source": "frozen ESM-2 esm2_t12_35M_UR50D mean-pooled embedding",
        }
    )
    result.to_csv(output_path, index=False)
    return result


def m0_shared_signal_correlations(
    wide: pd.DataFrame,
    output_path: Path,
    root: Path = ROOT,
) -> pd.DataFrame:
    sequence_path = root / "data" / "processed" / "v0_3_1" / "designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1.parquet"
    sequence = pd.read_parquet(sequence_path)[
        ["protein_id", "canonical_7mer", "kmer3_jaccard_to_paper_motif_rc_aware"]
    ].rename(
        columns={
            "protein_id": "dbp_id",
            "canonical_7mer": "canonical_7mer",
            "kmer3_jaccard_to_paper_motif_rc_aware": "kmer3_proxy",
        }
    )
    merged = wide.merge(sequence, on=["dbp_id", "canonical_7mer"], how="left", validate="one_to_one")
    merged["gc_fraction"] = merged["canonical_7mer"].map(lambda value: sum(base in "GC" for base in value) / len(value))
    rows = []
    for protein, group in merged.groupby("dbp_id", sort=True):
        rows.append(
            {
                "dbp_id": protein,
                "m0_vs_experimental": safe_spearman(group["M0_score"], group["experimental_score"]),
                "m0_vs_gc_fraction": safe_spearman(group["M0_score"], group["gc_fraction"]),
                "m0_vs_kmer3_proxy": safe_spearman(group["M0_score"], group["kmer3_proxy"]),
                "experimental_vs_gc_fraction": safe_spearman(group["experimental_score"], group["gc_fraction"]),
                "n_rc_classes": int(len(group)),
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)
    return result


def assay_reference_control_sensitivity(
    benchmark: pd.DataFrame,
    targets: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    rows = []
    target_frame = targets.set_index("dbp_id")
    for protein, group in benchmark.groupby("protein_id", sort=True):
        assay_target = target_frame.loc[protein, "experimental_assay_reference"]
        for name, scorer in (
            ("TargetHamming", target_hamming_control),
            ("TargetEdit", target_edit_control),
            ("TargetKmerOverlap", target_kmer_overlap_control),
        ):
            scores = [scorer(assay_target, candidate) for candidate in group["candidate_dna"]]
            rows.append(
                {
                    "dbp_id": protein,
                    "control": name,
                    "spearman": safe_spearman(group["experimental_score"], scores),
                    "assay_reference": assay_target,
                    "analysis_status": "SECONDARY_POST_HOC_DIAGNOSTIC",
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)
    return result


def replay_validation(
    replay_predictions: pd.DataFrame,
    primary_seed_results: pd.DataFrame,
    benchmark: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    rows = []
    for (fold_id, seed, protein, model), group in replay_predictions.groupby(
        ["fold_id", "seed", "dbp_id", "prediction_model"], sort=True
    ):
        truth = benchmark.loc[benchmark["protein_id"].eq(protein), "experimental_score"].to_numpy(dtype=float)
        prediction = group.sort_values("canonical_7mer")["prediction_score"].to_numpy(dtype=float)
        benchmark_group = benchmark.loc[benchmark["protein_id"].eq(protein)].sort_values("candidate_dna")
        metrics = compute_ranking_metrics(
            pd.DataFrame({"truth": benchmark_group["experimental_score"].to_numpy(), "prediction": prediction}),
            "truth",
            "prediction",
        )
        expected = primary_seed_results.loc[
            primary_seed_results["fold_id"].eq(fold_id)
            & primary_seed_results["seed"].eq(seed)
            & primary_seed_results["dbp_id"].eq(protein)
            & primary_seed_results["model"].eq(model)
        ]
        expected_value = float(expected["spearman"].iloc[0])
        rows.append(
            {
                "fold_id": fold_id,
                "seed": int(seed),
                "dbp_id": protein,
                "model": model,
                "replayed_spearman": metrics.spearman,
                "primary_spearman": expected_value,
                "absolute_difference": abs(metrics.spearman - expected_value),
                "matches_frozen_primary": bool(np.isclose(metrics.spearman, expected_value, atol=1e-10, rtol=1e-10)),
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(output_path, index=False)
    return result
