from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from src.v0_5_models import (
    CandidateDNAOnly,
    ProteinCandidate,
    ProteinTargetCandidate,
    TargetCandidateOnly,
    model_parameter_counts,
)
from src.v0_5_training import (
    V05Config,
    build_cache,
    evaluate_predictions,
    evaluate_target_controls,
    load_v05_data,
    predict_model,
    sample_rank_pairs,
    set_seed,
    train_model,
)


MODEL_ORDER = ("M0", "M1", "M1c", "M2", "M3")
MODEL_FULL_NAMES = {
    "M0": CandidateDNAOnly.model_name,
    "M1": ProteinCandidate.model_name,
    "M1c": "M1c_ProteinCandidateCapacityMatched",
    "M2": TargetCandidateOnly.model_name,
    "M3": ProteinTargetCandidate.model_name,
}
MODEL_INPUTS = {
    "M0": "D",
    "M1": "P,D",
    "M1c": "P,D",
    "M2": "T,D",
    "M3": "P,T,D",
}
DEVELOPMENT_EXPOSED = {"DBP1", "DBP3"}
UNSEEN_PRIMARY = {"DBP5", "DBP35", "DBP48", "DBP6", "DBP9"}


def model_factory(short_name: str, protein_dim: int, config: V05Config):
    if short_name == "M0":
        return CandidateDNAOnly(config.hidden_dim)
    if short_name == "M1":
        return ProteinCandidate(protein_dim, config.hidden_dim)
    if short_name == "M1c":
        return ProteinCandidate(protein_dim, config.capacity_matched_hidden_dim)
    if short_name == "M2":
        return TargetCandidateOnly(config.hidden_dim)
    if short_name == "M3":
        return ProteinTargetCandidate(protein_dim, config.hidden_dim)
    raise ValueError(f"Unknown v0.5 model: {short_name}")


def canonical_model_name(short_name: str) -> str:
    return ProteinCandidate.model_name if short_name == "M1c" else MODEL_FULL_NAMES[short_name]


def validate_loco_manifest(
    splits: pd.DataFrame,
    split_name: str,
    group_column: str,
) -> dict[str, Any]:
    """Validate that every DBP is held out once and leakage groups stay intact."""
    subset = splits.loc[splits["split_name"].eq(split_name)].copy()
    if subset.empty:
        raise ValueError(f"No rows found for split {split_name!r}")
    fold_ids = sorted(subset["fold_id"].unique())
    test_by_fold: dict[str, set[str]] = {}
    all_proteins = set(subset["dbp_id"])
    for fold_id in fold_ids:
        fold = subset.loc[subset["fold_id"].eq(fold_id)]
        if fold["dbp_id"].duplicated().any():
            raise AssertionError(f"Duplicate DBP rows in {split_name}/{fold_id}")
        train = set(fold.loc[fold["partition"].eq("train"), "dbp_id"])
        test = set(fold.loc[fold["partition"].eq("test"), "dbp_id"])
        if not train or not test or train & test:
            raise AssertionError(f"Invalid train/test partition in {split_name}/{fold_id}")
        for _, group in fold.groupby(group_column):
            if group["partition"].nunique() != 1:
                raise AssertionError(
                    f"{group_column} leakage in {split_name}/{fold_id}: "
                    f"{group[group_column].iloc[0]}"
                )
        test_by_fold[fold_id] = test
    coverage = set().union(*test_by_fold.values())
    counts = pd.Series([protein for test in test_by_fold.values() for protein in test]).value_counts()
    if coverage != all_proteins or not counts.eq(1).all():
        raise AssertionError(f"Test-fold coverage is not exactly once for {split_name}")
    return {
        "split_name": split_name,
        "fold_ids": fold_ids,
        "test_by_fold": test_by_fold,
        "all_proteins": sorted(all_proteins),
        "n_folds": len(fold_ids),
    }


def get_fold_partitions(
    splits: pd.DataFrame,
    split_name: str,
    fold_id: str,
) -> tuple[list[str], list[str]]:
    fold = splits.loc[splits["split_name"].eq(split_name) & splits["fold_id"].eq(fold_id)]
    if fold.empty:
        raise ValueError(f"Fold not found: {split_name}/{fold_id}")
    train = sorted(fold.loc[fold["partition"].eq("train"), "dbp_id"].unique())
    test = sorted(fold.loc[fold["partition"].eq("test"), "dbp_id"].unique())
    if set(train) & set(test):
        raise AssertionError(f"Protein leakage in {split_name}/{fold_id}")
    return train, test


def _seed_level_row(
    *,
    split_type: str,
    fold_id: str,
    seed: int,
    dbp_id: str,
    fold_metadata: pd.DataFrame,
    model: str,
    metrics: dict[str, Any],
    training_proteins: list[str],
    test_proteins: list[str],
    runtime: float,
    pair_count: int,
    status: str,
) -> dict[str, Any]:
    fold_row = fold_metadata.loc[fold_metadata["dbp_id"].eq(dbp_id)].iloc[0]
    return {
        "split_type": split_type,
        "fold_id": fold_id,
        "seed": int(seed),
        "dbp_id": dbp_id,
        "protein_cluster": fold_row["protein_cluster"],
        "combined_component": fold_row["combined_component"],
        "model": model,
        "n_units": int(metrics["n_rc_classes"]),
        "spearman": float(metrics["spearman"]),
        "ndcg_1pct": float(metrics["ndcg_1pct"]),
        "ndcg_5pct": float(metrics["ndcg_5pct"]),
        "pairwise_accuracy": float(metrics["pairwise_accuracy"]),
        "top1pct_recovery": float(metrics["top1pct_recovery"]),
        "training_proteins": "|".join(training_proteins),
        "test_proteins": "|".join(test_proteins),
        "runtime_seconds": float(runtime),
        "training_pair_count": int(pair_count),
        "development_exposed": dbp_id in DEVELOPMENT_EXPOSED,
        "status": status,
    }


def run_loco_evaluation(
    config: V05Config,
    *,
    split_name: str,
    group_column: str,
    seeds: Iterable[int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run all frozen models over a complete LOCO split and seed set."""
    benchmark, targets, embedding_map, splits = load_v05_data()
    validate_loco_manifest(splits, split_name, group_column)
    protein_dim = len(next(iter(embedding_map.values())))
    seed_rows: list[dict[str, Any]] = []
    health_rows: list[dict[str, Any]] = []
    for fold_id in sorted(splits.loc[splits["split_name"].eq(split_name), "fold_id"].unique()):
        train_proteins, test_proteins = get_fold_partitions(splits, split_name, fold_id)
        cache = build_cache(
            benchmark,
            targets,
            embedding_map,
            sorted(set(train_proteins + test_proteins)),
        )
        fold_metadata = splits.loc[
            splits["split_name"].eq(split_name) & splits["fold_id"].eq(fold_id)
        ]
        for seed in tuple(seeds):
            run_config = replace(config, seed=int(seed))
            pairs = sample_rank_pairs(
                benchmark,
                train_proteins,
                pairs_per_protein=config.pair_count_per_protein,
                seed=int(seed),
                tie_tolerance=config.tie_tolerance,
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
                predictions = {
                    canonical_model_name(short_name): {
                        protein: predict_model(model, canonical_model_name(short_name), cache[protein])
                        for protein in test_proteins
                    }
                }
                evaluated = evaluate_predictions(benchmark, predictions, test_proteins)
                for _, metrics in evaluated.iterrows():
                    seed_rows.append(
                        _seed_level_row(
                            split_type=split_name,
                            fold_id=fold_id,
                            seed=int(seed),
                            dbp_id=str(metrics["protein_id"]),
                            fold_metadata=fold_metadata,
                            model=short_name,
                            metrics=metrics.to_dict(),
                            training_proteins=train_proteins,
                            test_proteins=test_proteins,
                            runtime=runtime,
                            pair_count=len(pairs),
                            status="complete",
                        )
                    )
                all_predictions = np.concatenate(
                    [values for values in predictions[canonical_model_name(short_name)].values()]
                )
                first_loss = float(history["mean_pairwise_loss"].iloc[0])
                last_loss = float(history["mean_pairwise_loss"].iloc[-1])
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
                        "first_epoch_loss": first_loss,
                        "last_epoch_loss": last_loss,
                        "loss_delta": last_loss - first_loss,
                        "prediction_variance": float(np.var(all_predictions)),
                        "prediction_min": float(np.min(all_predictions)),
                        "prediction_max": float(np.max(all_predictions)),
                        "nan_inf_count": int(np.sum(~np.isfinite(all_predictions))),
                        "status": "complete",
                    }
                )
    return pd.DataFrame(seed_rows), pd.DataFrame(health_rows)


def aggregate_seed_level(seed_level: pd.DataFrame) -> pd.DataFrame:
    """Aggregate seeds within each protein/model before any macro summary."""
    group_columns = [
        "split_type",
        "fold_id",
        "dbp_id",
        "protein_cluster",
        "combined_component",
        "model",
        "development_exposed",
    ]
    metric_columns = [
        "spearman",
        "ndcg_1pct",
        "ndcg_5pct",
        "pairwise_accuracy",
        "top1pct_recovery",
    ]
    aggregations: dict[str, list[str]] = {column: ["mean", "std", "min", "max"] for column in metric_columns}
    result = seed_level.groupby(group_columns, sort=True, as_index=False).agg(aggregations)
    result.columns = [
        "_".join(part for part in column if part)
        if isinstance(column, tuple)
        else column
        for column in result.columns
    ]
    return result.rename(columns={"spearman_mean": "spearman"})


def build_per_protein_results(seed_level: pd.DataFrame) -> pd.DataFrame:
    aggregated = aggregate_seed_level(seed_level)
    index_columns = [
        "split_type",
        "fold_id",
        "dbp_id",
        "protein_cluster",
        "combined_component",
        "development_exposed",
    ]
    mean_table = aggregated.pivot_table(
        index=index_columns,
        columns="model",
        values="spearman",
        aggfunc="first",
    ).reset_index()
    sd_table = aggregated.pivot_table(
        index=index_columns,
        columns="model",
        values="spearman_std",
        aggfunc="first",
    ).reset_index()
    sd_table = sd_table.rename(columns={model: f"{model}_seed_sd" for model in MODEL_ORDER})
    result = mean_table.merge(sd_table, on=index_columns, how="left")
    for model in MODEL_ORDER:
        if model not in result.columns:
            result[model] = np.nan
        if f"{model}_seed_sd" not in result.columns:
            result[f"{model}_seed_sd"] = np.nan
    result["delta_m3_minus_m1"] = result["M3"] - result["M1"]
    result["delta_m3_minus_m1c"] = result["M3"] - result["M1c"]
    result["delta_m3_minus_m2"] = result["M3"] - result["M2"]
    ordered = index_columns + list(MODEL_ORDER) + [
        f"{model}_seed_sd" for model in MODEL_ORDER
    ] + [
        "delta_m3_minus_m1",
        "delta_m3_minus_m1c",
        "delta_m3_minus_m2",
    ]
    return result[ordered].sort_values(["split_type", "fold_id", "dbp_id"]).reset_index(drop=True)


def _macro_stats(values: pd.Series) -> tuple[float, float, float]:
    clean = values.dropna().to_numpy(dtype=float)
    if not len(clean):
        return np.nan, np.nan, np.nan
    return float(np.median(clean)), float(np.mean(clean)), float(np.std(clean, ddof=1)) if len(clean) > 1 else np.nan


def build_primary_macro_summary(per_protein: pd.DataFrame, seeds: Iterable[int]) -> pd.DataFrame:
    all_rows = per_protein.loc[per_protein["dbp_id"].isin(sorted(set(per_protein["dbp_id"])))]
    unseen_rows = per_protein.loc[per_protein["dbp_id"].isin(UNSEEN_PRIMARY)]
    delta_columns = {
        "median_delta_m3_minus_m1": "delta_m3_minus_m1",
        "median_delta_m3_minus_m1c": "delta_m3_minus_m1c",
        "median_delta_m3_minus_m2": "delta_m3_minus_m2",
    }
    improvement_columns = {
        "m3_improved_over_m1": "delta_m3_minus_m1",
        "m3_improved_over_m1c": "delta_m3_minus_m1c",
        "m3_improved_over_m2": "delta_m3_minus_m2",
    }
    rows = []
    for model in MODEL_ORDER:
        all_median, all_mean, all_sd = _macro_stats(all_rows[model])
        unseen_median, unseen_mean, unseen_sd = _macro_stats(unseen_rows[model])
        row: dict[str, Any] = {
            "model": model,
            "all7_macro_median": all_median,
            "unseen5_macro_median": unseen_median,
            "all7_macro_mean": all_mean,
            "unseen5_macro_mean": unseen_mean,
            "all7_macro_sd": all_sd,
            "unseen5_macro_sd": unseen_sd,
            "proteins_evaluated": int(all_rows[model].notna().sum()),
            "unseen5_proteins_evaluated": int(unseen_rows[model].notna().sum()),
            "seeds": "|".join(str(seed) for seed in seeds),
        }
        for output_name, source_column in delta_columns.items():
            all_value, _, _ = _macro_stats(all_rows[source_column])
            unseen_value, _, _ = _macro_stats(unseen_rows[source_column])
            row[f"{output_name}_all7"] = all_value
            row[f"{output_name}_unseen5"] = unseen_value
        for output_name, source_column in improvement_columns.items():
            row[f"{output_name}_all7_count"] = int((all_rows[source_column] > 0).sum())
            row[f"{output_name}_unseen5_count"] = int((unseen_rows[source_column] > 0).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def build_strict_macro_summary(per_protein: pd.DataFrame, seeds: Iterable[int]) -> pd.DataFrame:
    rows = []
    for model in MODEL_ORDER:
        median, mean, sd = _macro_stats(per_protein[model])
        rows.append(
            {
                "model": model,
                "all_evaluated_macro_median": median,
                "all_evaluated_macro_mean": mean,
                "all_evaluated_macro_sd": sd,
                "proteins_evaluated": int(per_protein[model].notna().sum()),
                "components_evaluated": int(per_protein["fold_id"].nunique()),
                "seeds": "|".join(str(seed) for seed in seeds),
                "median_delta_m3_minus_m1": _macro_stats(per_protein["delta_m3_minus_m1"])[0],
                "median_delta_m3_minus_m1c": _macro_stats(per_protein["delta_m3_minus_m1c"])[0],
                "median_delta_m3_minus_m2": _macro_stats(per_protein["delta_m3_minus_m2"])[0],
                "m3_improved_over_m1_count": int((per_protein["delta_m3_minus_m1"] > 0).sum()),
                "m3_improved_over_m1c_count": int((per_protein["delta_m3_minus_m1c"] > 0).sum()),
                "m3_improved_over_m2_count": int((per_protein["delta_m3_minus_m2"] > 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def build_baseline_context(
    primary_per_protein: pd.DataFrame,
    *,
    root: Path,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    def add(
        baseline: str,
        coverage: str,
        spearman: float,
        regime: str,
        matched: bool,
        notes: str,
    ) -> None:
        rows.append(
            {
                "baseline": baseline,
                "coverage": coverage,
                "macro_median_spearman": spearman,
                "split_training_regime": regime,
                "directly_matched": matched,
                "notes": notes,
            }
        )

    sequence = pd.read_csv(
        root / "results" / "v0_3_1" / "tables" / "designed_dbp_sequence_baseline_rc_aware.csv"
    )
    sequence_values = sequence.loc[sequence["metric"].eq("kmer3_jaccard_to_paper_motif_rc_aware"), "spearman"]
    add(
        "sequence-only kmer3 RC-aware",
        "7/7",
        float(sequence_values.median()),
        "v0.3.1 designed uPBM; no training",
        True,
        "Sequence-only proxy; not protein-conditioned.",
    )

    v041 = pd.read_csv(root / "results" / "v0_4_1" / "tables" / "v0_4_1_baseline_summary.csv")
    simple = v041.loc[
        v041["baseline"].str.contains("SimpleProteinConditionalBaseline")
        & v041["evaluation_dataset"].str.contains("designed"),
        "macro_median_spearman",
    ]
    natural_simple = v041.loc[
        v041["baseline"].str.contains("SimpleProteinConditionalBaseline")
        & v041["evaluation_dataset"].str.contains("natural_test"),
        "macro_median_spearman",
    ]
    add(
        "SimpleProteinConditional",
        "7/7 designed; 9 natural test",
        float(simple.iloc[0]),
        "natural PBM train -> designed external; prior v0.4.1",
        False,
        f"Prior natural-test macro median={float(natural_simple.iloc[0]):.6f}; low-capacity baseline.",
    )

    v042 = pd.read_csv(
        root / "results" / "v0_4_2" / "tables" / "final_strong_baseline_summary_deeppbs_completed_v0_4_2.csv"
    )
    frozen = v042.loc[v042["method"].str.contains("FrozenPLM"), "designed_macro_median_spearman"]
    add(
        "FrozenPLM",
        "7/7 designed",
        float(frozen.iloc[0]),
        "natural PBM train -> designed external; prior v0.4.2",
        False,
        "Frozen ESM-2 baseline; prior result is not this matched v0.5 trainer.",
    )
    deeppbs = v042.loc[v042["method"].eq("DeepPBS"), "designed_macro_median_spearman"]
    add(
        "DeepPBS",
        "2/7 designed",
        float(deeppbs.iloc[0]),
        "official structure-aware diagnostic; prior v0.4.2",
        False,
        "Coverage-limited diagnostic; not comparable as a seven-protein estimate.",
    )
    nampnn = v042.loc[v042["method"].str.contains("NA-MPNN"), "designed_macro_median_spearman"]
    add(
        "NA-MPNN diagnostic",
        "2/7 designed",
        float(nampnn.iloc[0]),
        "official structure-aware diagnostic; prior v0.4.2",
        False,
        "Limited diagnostic coverage and DBP48 overlap caveat.",
    )
    replicate = v042.loc[v042["method"].str.contains("Replicate"), "designed_macro_median_spearman"]
    add(
        "Experimental replicate reference",
        "7/7 designed",
        float(replicate.iloc[0]),
        "uPBM replicate agreement",
        False,
        "Empirical reproducibility reference, not a strict theoretical ceiling.",
    )

    for model in MODEL_ORDER:
        value = float(primary_per_protein[model].median())
        add(
            f"v0.5 {model}",
            "7/7 designed",
            value,
            "v0.5 4-fold protein-cluster LOCO, seed-mean per protein",
            True,
            "Matched primary result; see primary tables for development-exposed and unseen subsets.",
        )
    return pd.DataFrame(rows)


def build_primary_artifacts(config: V05Config, root: Path) -> dict[str, pd.DataFrame]:
    seeds = tuple(config.evaluation_seeds)
    primary_seed, primary_health = run_loco_evaluation(
        config,
        split_name="protein_cluster_loco",
        group_column="protein_cluster",
        seeds=seeds,
    )
    strict_seed, strict_health = run_loco_evaluation(
        config,
        split_name="combined_component_loco",
        group_column="combined_component",
        seeds=seeds,
    )
    primary_per_protein = build_per_protein_results(primary_seed)
    strict_per_protein = build_per_protein_results(strict_seed)
    primary_macro = build_primary_macro_summary(primary_per_protein, seeds)
    strict_macro = build_strict_macro_summary(strict_per_protein, seeds)
    benchmark, targets, _, _ = load_v05_data()
    controls = evaluate_target_controls(
        benchmark,
        targets,
        sorted(benchmark["protein_id"].unique()),
    )
    controls.insert(0, "split_scope", "all_7_proteins")
    context = build_baseline_context(primary_per_protein, root=root)
    health = pd.concat([primary_health, strict_health], ignore_index=True)
    return {
        "primary_seed": primary_seed,
        "primary_health": primary_health,
        "strict_seed": strict_seed,
        "strict_health": strict_health,
        "primary_per_protein": primary_per_protein,
        "strict_per_protein": strict_per_protein,
        "primary_macro": primary_macro,
        "strict_macro": strict_macro,
        "controls": controls,
        "baseline_context": context,
        "training_health": health,
    }
