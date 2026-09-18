from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[2]
SAMBA_SOURCE = (
    ROOT
    / "data/raw/v1_3_2_external/samba/41586_2020_2843_MOESM7_ESM.xlsx"
)
SPOT_TABLE = ROOT / "data/processed/v1_3_2_samba_spots.parquet"
RESULTS = ROOT / "results/v1_3"
FIGURES = ROOT / "figures/v1_3"
SEED = 1301
N_PERM = 2_000

SITE_SHEETS = [
    "Ets1_site_1",
    "Ets1_site_2",
    "Ets1_site_3",
    "Ets1_site_4",
    "TBP_site_1",
    "TBP_site_2",
    "Max_site_1",
    "Max_site_2",
    "Cbf1",
    "Egr1",
    "p53",
    "GR",
]

RAW_LAYOUT = {
    "Raw_data1": {
        "p53": "p53_Signal",
        "Ets1": "Ets1_Signal",
        "Cbf1": "Cbf1_Signal",
        "GR": "GR_Signal",
    },
    "Raw_data2": {"Max": "Max_Signal", "Egr1": "Egr1_Signal"},
    "Raw_data3": {"TBP": "TBP_Signal"},
}


def safe_corr(x: np.ndarray, y: np.ndarray, kind: str = "spearman") -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return math.nan
    fn = spearmanr if kind == "spearman" else pearsonr
    return float(fn(x, y).statistic)


def pairwise_accuracy(experimental: np.ndarray, estimate: np.ndarray) -> tuple[float, int]:
    experimental = np.asarray(experimental, dtype=float)
    estimate = np.asarray(estimate, dtype=float)
    scores: list[float] = []
    for left in range(len(experimental)):
        for right in range(left + 1, len(experimental)):
            exp_diff = experimental[left] - experimental[right]
            est_diff = estimate[left] - estimate[right]
            if exp_diff == 0:
                continue
            scores.append(1.0 if exp_diff * est_diff > 0 else 0.5 if est_diff == 0 else 0.0)
    return (float(np.mean(scores)), len(scores)) if scores else (math.nan, 0)


def tf_for_site(site: str) -> str:
    return site.split("_", 1)[0]


def wt_basepair(site_table: pd.DataFrame, position: int) -> str:
    wt = site_table.loc[site_table["Position"].astype(str).eq("WT")].iloc[0]
    strand1 = str(wt["Strand1"])
    strand2 = str(wt["Strand2"])
    if not (1 <= position <= len(strand1) == len(strand2)):
        raise ValueError(f"Invalid position {position} for {strand1}/{strand2}")
    return strand1[position - 1] + strand2[-position]


def load_source_tables() -> tuple[dict[str, pd.DataFrame], dict[str, tuple[str, pd.DataFrame, str]]]:
    requested = SITE_SHEETS + list(RAW_LAYOUT)
    tables = pd.read_excel(SAMBA_SOURCE, sheet_name=requested)
    sites = {name: tables[name] for name in SITE_SHEETS}
    raw_by_tf: dict[str, tuple[str, pd.DataFrame, str]] = {}
    for sheet, tf_columns in RAW_LAYOUT.items():
        raw = tables[sheet]
        for tf_id, signal_column in tf_columns.items():
            raw_by_tf[tf_id] = (sheet, raw, signal_column)
    return sites, raw_by_tf


def recover_spot_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not SAMBA_SOURCE.exists():
        raise FileNotFoundError(SAMBA_SOURCE)
    sites, raw_by_tf = load_source_tables()
    rows: list[dict] = []
    audit_rows: list[dict] = []

    for site_name in SITE_SHEETS:
        site = sites[site_name].copy()
        tf_id = tf_for_site(site_name)
        raw_sheet, raw, signal_column = raw_by_tf[tf_id]
        raw_subset = raw[["ID", "Name", "Sequence", signal_column]].copy()
        raw_subset[signal_column] = pd.to_numeric(raw_subset[signal_column], errors="coerce")
        raw_subset = raw_subset.dropna(subset=["Sequence", signal_column])

        mapping: list[dict] = []
        for source_row, variant in site.iterrows():
            is_reference = str(variant["Position"]) == "WT"
            position = None if is_reference else int(variant["Position"])
            perturbation_type = (
                "reference"
                if is_reference
                else "mismatch"
                if int(variant["Is_mismatch"]) == 1
                else "canonical_watson_crick"
            )
            mapping.append(
                {
                    "sequence": str(variant["Sequence"]),
                    "position": position,
                    "wt_basepair": "WT" if is_reference else wt_basepair(site, position),
                    "perturbation": "WT" if is_reference else str(variant["String"]),
                    "perturbation_type": perturbation_type,
                    "published_median_signal": float(variant["Median_signal"]),
                    "published_fold_change": float(variant["Fold_change"]),
                    "published_calibrated_dg": float(variant["Calibrated_dG"]),
                    "source_site_row": int(source_row + 2),
                }
            )
        mapping_table = pd.DataFrame(mapping)
        if mapping_table.duplicated(["sequence"]).any():
            duplicates = mapping_table.loc[mapping_table.duplicated(["sequence"], keep=False), "sequence"]
            raise ValueError(f"Duplicate mapped sequence within {site_name}: {duplicates.iloc[0]}")

        joined = mapping_table.merge(
            raw_subset,
            left_on="sequence",
            right_on="Sequence",
            how="left",
            validate="one_to_many",
            indicator=True,
        )
        matched = joined[joined["_merge"].eq("both")].copy()
        missing = joined[joined["_merge"].eq("left_only")]

        # Some calibration libraries reuse an identical 60-mer in more than
        # one named block.  Published medians were sometimes calculated over
        # all such blocks, so a name collision alone is not an exclusion.
        # Validate the sequence join by reconstructing every published median;
        # never select a block by its agreement with an individual outcome.
        matched["name_prefix"] = matched["Name"].astype(str).str.replace(
            r"_r\d+$", "", regex=True
        )
        prefix_counts = matched.groupby("sequence")["name_prefix"].nunique()
        name_collision_sequences = set(prefix_counts[prefix_counts > 1].index)
        reconstructed = matched.groupby("sequence").agg(
            reconstructed_median=(signal_column, "median"),
            published_median=("published_median_signal", "first"),
        )
        reconstructed["absolute_error"] = (
            reconstructed.reconstructed_median - reconstructed.published_median
        ).abs()
        reconstruction_mismatches = set(
            reconstructed.index[reconstructed.absolute_error > 1e-12]
        )
        site_mapping_status = (
            "UNRESOLVED_PUBLISHED_MEDIAN_RECONSTRUCTION"
            if reconstruction_mismatches
            else "VALIDATED_BY_EXACT_MEDIAN_RECONSTRUCTION"
        )

        reference = matched[matched["perturbation_type"].eq("reference")][signal_column]
        reference_mean = float(reference.mean()) if len(reference) else math.nan
        if not np.isfinite(reference_mean) or reference_mean <= 0:
            raise ValueError(f"No positive reference signal for {site_name}/{tf_id}")

        for _, spot in matched.iterrows():
            raw_signal = float(spot[signal_column])
            rows.append(
                {
                    "tf_id": tf_id,
                    "binding_site_id": site_name,
                    "position": spot["position"],
                    "wt_basepair": spot["wt_basepair"],
                    "perturbation": spot["perturbation"],
                    "perturbation_type": spot["perturbation_type"],
                    "sequence": spot["sequence"],
                    "spot_id": f"{raw_sheet}:{tf_id}:{spot['ID']}",
                    "source_spot_name": str(spot["Name"]),
                    # The workbook exposes raw rows but no trustworthy array
                    # or batch identifier. Keep this explicit rather than
                    # treating the source sheet name as an experimental array.
                    "array_id": "unknown_array",
                    "array_id_source": "not_supplied_by_source_workbook",
                    "raw_signal": raw_signal,
                    "normalized_signal": math.log2(raw_signal / reference_mean),
                    "normalization_scope": "full_site_reference_mean_descriptive_only",
                    "published_median_signal": float(spot["published_median_signal"]),
                    "published_fold_change": float(spot["published_fold_change"]),
                    "published_calibrated_dg": float(spot["published_calibrated_dg"]),
                    "source_file": SAMBA_SOURCE.relative_to(ROOT).as_posix(),
                    "source_sheet": raw_sheet,
                    "source_site_sheet": site_name,
                    "source_site_row": int(spot["source_site_row"]),
                    "source_mapping_status": site_mapping_status,
                }
            )

        spot_counts = matched.groupby("perturbation_type").size().to_dict()
        sequence_counts = matched.groupby("perturbation_type")["sequence"].nunique().to_dict()
        audit_rows.append(
            {
                "tf_id": tf_id,
                "binding_site_id": site_name,
                "n_positions": int(mapping_table["position"].nunique()),
                "n_canonical_sequences": int(sequence_counts.get("canonical_watson_crick", 0)),
                "n_mismatch_sequences": int(sequence_counts.get("mismatch", 0)),
                "n_reference_sequences": int(sequence_counts.get("reference", 0)),
                "n_canonical_spots": int(spot_counts.get("canonical_watson_crick", 0)),
                "n_mismatch_spots": int(spot_counts.get("mismatch", 0)),
                "n_reference_spots": int(spot_counts.get("reference", 0)),
                "n_unmapped_sequences": int(missing["sequence"].nunique()),
                "n_name_collision_sequences": len(name_collision_sequences),
                "n_reconstruction_mismatches": len(reconstruction_mismatches),
                "max_published_median_reconstruction_error": float(reconstructed.absolute_error.max()),
                "median_published_median_reconstruction_error": float(reconstructed.absolute_error.median()),
                "source_mapping_status": site_mapping_status,
                "min_spots_per_mapped_sequence": int(matched.groupby("sequence").size().min()),
                "median_spots_per_mapped_sequence": float(matched.groupby("sequence").size().median()),
                "max_spots_per_mapped_sequence": int(matched.groupby("sequence").size().max()),
                "raw_sheet": raw_sheet,
                "signal_column": signal_column,
            }
        )

    spots = pd.DataFrame(rows).sort_values(
        ["tf_id", "binding_site_id", "perturbation_type", "position", "perturbation", "spot_id"],
        na_position="first",
    ).reset_index(drop=True)
    spots["position"] = spots["position"].astype("Int64")
    if spots.duplicated(["tf_id", "binding_site_id", "sequence", "spot_id"]).any():
        raise ValueError("Duplicate SaMBA spot key")
    audit = pd.DataFrame(audit_rows).sort_values(["tf_id", "binding_site_id"]).reset_index(drop=True)
    SPOT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    spots.to_parquet(SPOT_TABLE, index=False)
    audit.to_csv(RESULTS / "external_samba_data_audit.tsv", sep="\t", index=False)
    return spots, audit


def assign_split(spots: pd.DataFrame, seed: int) -> pd.DataFrame:
    result = spots.sort_values(
        ["tf_id", "binding_site_id", "sequence", "array_id", "spot_id"]
    ).reset_index(drop=True).copy()
    result["technical_split"] = ""
    rng = np.random.Generator(np.random.PCG64(seed))
    keys = ["tf_id", "binding_site_id", "sequence", "array_id"]
    for _, group in result.groupby(keys, sort=True, dropna=False):
        if len(group) < 2:
            continue
        order = rng.permutation(group.index.to_numpy())
        labels = np.where(np.arange(len(order)) % 2 == 0, "A", "B")
        result.loc[order, "technical_split"] = labels
    return result


def split_effects(spots: pd.DataFrame, seed: int) -> pd.DataFrame:
    assigned = assign_split(spots, seed)
    valid = assigned[
        assigned["technical_split"].isin(["A", "B"])
        & np.isfinite(assigned["raw_signal"])
        & assigned["raw_signal"].gt(0)
        & assigned["source_mapping_status"].eq("VALIDATED_BY_EXACT_MEDIAN_RECONSTRUCTION")
    ].copy()
    group_keys = [
        "tf_id",
        "binding_site_id",
        "position",
        "wt_basepair",
        "perturbation",
        "perturbation_type",
        "sequence",
        "array_id",
        "technical_split",
    ]
    estimates = (
        valid.groupby(group_keys, observed=True, dropna=False)["raw_signal"]
        .agg(mean_signal="mean", n_spots="size")
        .reset_index()
    )
    references = estimates[estimates["perturbation_type"].eq("reference")][
        ["tf_id", "binding_site_id", "array_id", "technical_split", "mean_signal", "n_spots"]
    ].rename(columns={"mean_signal": "reference_mean_signal", "n_spots": "reference_n_spots"})
    mutations = estimates[~estimates["perturbation_type"].eq("reference")].copy()
    mutations = mutations.merge(
        references,
        on=["tf_id", "binding_site_id", "array_id", "technical_split"],
        how="left",
        validate="many_to_one",
    )
    mutations["effect"] = np.log2(mutations["mean_signal"] / mutations["reference_mean_signal"])
    wide = mutations.pivot(
        index=[
            "tf_id",
            "binding_site_id",
            "position",
            "wt_basepair",
            "perturbation",
            "perturbation_type",
            "sequence",
            "array_id",
        ],
        columns="technical_split",
        values=["effect", "mean_signal", "n_spots", "reference_mean_signal", "reference_n_spots"],
    )
    wide.columns = [f"{value}_{split.lower()}" for value, split in wide.columns]
    wide = wide.reset_index()
    wide["split_seed"] = seed
    required = ["effect_a", "effect_b", "mean_signal_a", "mean_signal_b"]
    return wide.dropna(subset=required).sort_values(
        ["tf_id", "binding_site_id", "perturbation_type", "position", "perturbation"]
    ).reset_index(drop=True)


def landscape_variance_fraction(table: pd.DataFrame, value: str) -> tuple[float, float]:
    observed = table[value].to_numpy(float)
    if len(observed) == 0 or np.var(observed, ddof=0) == 0:
        return math.nan, math.nan
    position_mean = table.groupby("position")[value].transform("mean").to_numpy(float)
    alpha = position_mean - np.mean(observed)
    residual = observed - position_mean
    total = np.var(observed, ddof=0)
    return float(np.var(alpha, ddof=0) / total), float(np.var(residual, ddof=0) / total)


def evaluate_site(table: pd.DataFrame) -> dict:
    table = table.sort_values(["position", "perturbation"]).reset_index(drop=True)
    complete = table[np.isfinite(table["effect_a"]) & np.isfinite(table["effect_b"])].copy()
    result: dict[str, float | int | str] = {
        "n_perturbations": len(complete),
        "n_positions": int(complete["position"].nunique()),
        "global_status": "EVALUABLE" if len(complete) >= 6 else "NOT_EVALUABLE",
    }
    if len(complete) >= 6:
        result["global_spearman"] = safe_corr(complete.effect_a, complete.effect_b)
        result["global_pearson"] = safe_corr(complete.effect_a, complete.effect_b, "pearson")
        result["global_mae"] = float(np.mean(np.abs(complete.effect_a - complete.effect_b)))
    else:
        result.update(global_spearman=math.nan, global_pearson=math.nan, global_mae=math.nan)

    position = complete.groupby("position").agg(
        sensitivity_a=("effect_a", lambda values: float(np.mean(np.abs(values)))),
        sensitivity_b=("effect_b", lambda values: float(np.mean(np.abs(values)))),
    )
    position_ok = len(position) >= 3
    result["position_status"] = "EVALUABLE" if position_ok else "NOT_EVALUABLE"
    result["position_spearman"] = safe_corr(position.sensitivity_a, position.sensitivity_b) if position_ok else math.nan
    result["position_pearson"] = safe_corr(position.sensitivity_a, position.sensitivity_b, "pearson") if position_ok else math.nan

    eligible_positions = [
        position_id
        for position_id, group in complete.groupby("position", sort=True)
        if len(group) >= 2 and group["perturbation"].nunique() >= 2
    ]
    identity = complete[complete["position"].isin(eligible_positions)].copy()
    identity_ok = len(identity) >= 6 and len(eligible_positions) >= 3
    identity["residual_a"] = identity.effect_a - identity.groupby("position").effect_a.transform("mean")
    identity["residual_b"] = identity.effect_b - identity.groupby("position").effect_b.transform("mean")
    pair_scores: list[float] = []
    n_pairs = 0
    for _, group in identity.groupby("position", sort=True):
        score, count = pairwise_accuracy(group.effect_a.to_numpy(), group.effect_b.to_numpy())
        if count:
            pair_scores.extend([score] * count)
            n_pairs += count
    pair_ok = identity_ok and n_pairs >= 3
    result["identity_status"] = "EVALUABLE" if identity_ok else "NOT_EVALUABLE"
    result["identity_pairwise_status"] = "EVALUABLE" if pair_ok else "NOT_EVALUABLE"
    result["n_identity_perturbations"] = len(identity)
    result["n_identity_positions"] = len(eligible_positions)
    result["n_identity_pairs"] = n_pairs
    result["identity_residual_spearman"] = safe_corr(identity.residual_a, identity.residual_b) if identity_ok else math.nan
    result["identity_residual_pearson"] = safe_corr(identity.residual_a, identity.residual_b, "pearson") if identity_ok else math.nan
    if pair_ok:
        # Recompute from all pairs so positions with more comparable pairs receive their natural row weight.
        correct = []
        for _, group in identity.groupby("position", sort=True):
            values = group[["effect_a", "effect_b"]].to_numpy(float)
            for left in range(len(values)):
                for right in range(left + 1, len(values)):
                    exp_diff = values[left, 0] - values[right, 0]
                    est_diff = values[left, 1] - values[right, 1]
                    if exp_diff == 0:
                        continue
                    correct.append(1.0 if exp_diff * est_diff > 0 else 0.5 if est_diff == 0 else 0.0)
        pair_accuracy = float(np.mean(correct))
    else:
        pair_accuracy = math.nan
    result["identity_pairwise_accuracy"] = pair_accuracy
    result["identity_kendall_style_concordance"] = 2 * pair_accuracy - 1 if np.isfinite(pair_accuracy) else math.nan

    if identity_ok:
        pos_a, id_a = landscape_variance_fraction(identity, "effect_a")
        pos_b, id_b = landscape_variance_fraction(identity, "effect_b")
    else:
        pos_a = id_a = pos_b = id_b = math.nan
    result.update(
        fraction_position_a=pos_a,
        fraction_identity_a=id_a,
        fraction_position_b=pos_b,
        fraction_identity_b=id_b,
    )
    return result


def primary_decomposition(effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, group in effects.groupby(
        ["tf_id", "binding_site_id", "perturbation_type"], sort=True
    ):
        tf_id, site_id, perturbation_type = keys
        rows.append(
            {
                "tf_id": tf_id,
                "binding_site_id": site_id,
                "perturbation_type": perturbation_type,
                "comparison": "technical_spot_split_A_vs_B",
                "split_seed": SEED,
                **evaluate_site(group),
            }
        )
    result = pd.DataFrame(rows).sort_values(
        ["perturbation_type", "tf_id", "binding_site_id"]
    ).reset_index(drop=True)
    result.to_csv(RESULTS / "external_samba_decomposition.tsv", sep="\t", index=False)
    result[result.perturbation_type.eq("canonical_watson_crick")].to_csv(
        RESULTS / "external_samba_canonical_decomposition.tsv", sep="\t", index=False
    )
    result[result.perturbation_type.eq("mismatch")].to_csv(
        RESULTS / "external_samba_mismatch_decomposition.tsv", sep="\t", index=False
    )
    return result


def secondary_split_stability(spots: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for seed in range(SEED, SEED + 100):
        effects = split_effects(spots, seed)
        for keys, group in effects.groupby(
            ["tf_id", "binding_site_id", "perturbation_type"], sort=True
        ):
            tf_id, site_id, perturbation_type = keys
            metrics = evaluate_site(group)
            rows.append(
                {
                    "tf_id": tf_id,
                    "binding_site_id": site_id,
                    "perturbation_type": perturbation_type,
                    "split_seed": seed,
                    "global_spearman": metrics["global_spearman"],
                    "position_spearman": metrics["position_spearman"],
                    "identity_residual_spearman": metrics["identity_residual_spearman"],
                    "identity_pairwise_accuracy": metrics["identity_pairwise_accuracy"],
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(RESULTS / "external_samba_split_stability.tsv", sep="\t", index=False)
    return result


def permutation_nulls(effects: pd.DataFrame, observed: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    distribution_rows = []
    rng = np.random.default_rng(SEED)
    observed_index = observed.set_index(["tf_id", "binding_site_id", "perturbation_type"])
    for keys, group in effects.groupby(
        ["tf_id", "binding_site_id", "perturbation_type"], sort=True
    ):
        tf_id, site_id, perturbation_type = keys
        group = group.sort_values(["position", "perturbation"]).reset_index(drop=True)
        obs = observed_index.loc[keys]
        positions = [part.index.to_numpy() for _, part in group.groupby("position", sort=True)]
        pos_a = group.groupby("position").effect_a.apply(lambda x: float(np.mean(np.abs(x)))).to_numpy()
        pos_b = group.groupby("position").effect_b.apply(lambda x: float(np.mean(np.abs(x)))).to_numpy()
        values = {
            "position_spearman": [],
            "identity_residual_spearman": [],
            "identity_pairwise_accuracy": [],
        }
        for _ in range(N_PERM):
            values["position_spearman"].append(safe_corr(pos_a, rng.permutation(pos_b)))
            shuffled = group.effect_b.to_numpy(float).copy()
            for indices in positions:
                shuffled[indices] = rng.permutation(shuffled[indices])
            residual_a = group.effect_a - group.groupby("position").effect_a.transform("mean")
            shuffled_series = pd.Series(shuffled, index=group.index)
            residual_b = shuffled_series - shuffled_series.groupby(group.position).transform("mean")
            values["identity_residual_spearman"].append(safe_corr(residual_a, residual_b))
            pair_values = []
            for indices in positions:
                score, count = pairwise_accuracy(group.effect_a.iloc[indices], shuffled[indices])
                if count:
                    pair_values.extend([score] * count)
            values["identity_pairwise_accuracy"].append(
                float(np.mean(pair_values)) if pair_values else math.nan
            )

        definitions = [
            ("position_label", "position_spearman", float(obs.position_spearman), 0.0, "two_sided"),
            (
                "within_position_base_label",
                "identity_residual_spearman",
                float(obs.identity_residual_spearman),
                0.0,
                "two_sided",
            ),
            (
                "within_position_base_label",
                "identity_pairwise_accuracy",
                float(obs.identity_pairwise_accuracy),
                0.5,
                "upper",
            ),
        ]
        for null_type, metric, observed_value, chance, tail in definitions:
            null = np.asarray(values[metric], dtype=float)
            null = null[np.isfinite(null)]
            if not np.isfinite(observed_value) or len(null) == 0:
                empirical_p = math.nan
            elif tail == "two_sided":
                empirical_p = float(
                    (1 + np.sum(np.abs(null - chance) >= abs(observed_value - chance)))
                    / (1 + len(null))
                )
            else:
                empirical_p = float((1 + np.sum(null >= observed_value)) / (1 + len(null)))
            summary_rows.append(
                {
                    "tf_id": tf_id,
                    "binding_site_id": site_id,
                    "perturbation_type": perturbation_type,
                    "null_type": null_type,
                    "metric": metric,
                    "observed": observed_value,
                    "chance_baseline": chance,
                    "null_median": float(np.median(null)) if len(null) else math.nan,
                    "null_q025": float(np.quantile(null, 0.025)) if len(null) else math.nan,
                    "null_q975": float(np.quantile(null, 0.975)) if len(null) else math.nan,
                    "empirical_p": empirical_p,
                    "tail": tail,
                    "n_permutations": len(null),
                }
            )
            distribution_rows.extend(
                {
                    "tf_id": tf_id,
                    "binding_site_id": site_id,
                    "perturbation_type": perturbation_type,
                    "null_type": null_type,
                    "metric": metric,
                    "permutation": permutation,
                    "null_value": value,
                }
                for permutation, value in enumerate(null)
            )
    summary = pd.DataFrame(summary_rows)
    distribution = pd.DataFrame(distribution_rows)
    summary.to_csv(RESULTS / "external_samba_permutation_nulls.tsv", sep="\t", index=False)
    distribution.to_parquet(RESULTS / "external_samba_permutation_distributions.parquet", index=False)
    return summary, distribution


def tf_level_summary(decomposition: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics = [
        "global_spearman",
        "position_spearman",
        "identity_residual_spearman",
        "identity_pairwise_accuracy",
        "fraction_position_a",
        "fraction_identity_a",
        "fraction_position_b",
        "fraction_identity_b",
    ]
    per_tf = (
        decomposition.groupby(["perturbation_type", "tf_id"], sort=True)[metrics]
        .median()
        .reset_index()
    )
    rows = []
    rng = np.random.default_rng(SEED)
    for perturbation_type, group in per_tf.groupby("perturbation_type", sort=True):
        for metric in metrics:
            values = group[metric].dropna().to_numpy(float)
            if len(values) >= 5:
                boot = np.asarray(
                    [np.median(rng.choice(values, len(values), replace=True)) for _ in range(2_000)]
                )
                low, high = np.quantile(boot, [0.025, 0.975])
            else:
                low = high = math.nan
            rows.append(
                {
                    "perturbation_type": perturbation_type,
                    "metric": metric,
                    "n_tf": len(values),
                    "median": float(np.median(values)) if len(values) else math.nan,
                    "mean": float(np.mean(values)) if len(values) else math.nan,
                    "ci_low": low,
                    "ci_high": high,
                    "aggregation": "site_median_within_tf_then_across_tf",
                }
            )
    summary = pd.DataFrame(rows)
    per_tf.to_csv(RESULTS / "external_samba_decomposition_per_tf.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "external_samba_decomposition_summary.tsv", sep="\t", index=False)
    return per_tf, summary


def render_figure3(per_tf: pd.DataFrame) -> pd.DataFrame:
    metric_labels = {
        "global_spearman": "Global",
        "position_spearman": "Position",
        "identity_residual_spearman": "Identity residual",
        "identity_pairwise_accuracy": "Identity pairwise",
    }
    source = per_tf.melt(
        id_vars=["perturbation_type", "tf_id"],
        value_vars=list(metric_labels),
        var_name="metric",
        value_name="value",
    )
    source["metric_label"] = source.metric.map(metric_labels)
    source["chance_baseline"] = np.where(source.metric.eq("identity_pairwise_accuracy"), 0.5, 0.0)
    source["analysis_role"] = np.where(
        source.perturbation_type.eq("canonical_watson_crick"), "PRIMARY", "SECONDARY"
    )
    source.to_csv(RESULTS / "figure3_external_identifiability.tsv", sep="\t", index=False)

    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper")
    order = list(metric_labels.values())
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1), sharey=True)
    panels = [
        ("canonical_watson_crick", "Canonical Watson-Crick mutations (primary)"),
        ("mismatch", "Mismatches (secondary)"),
    ]
    for axis, (perturbation_type, title) in zip(axes, panels, strict=True):
        panel = source[source.perturbation_type.eq(perturbation_type)].copy()
        sns.stripplot(
            data=panel,
            x="metric_label",
            y="value",
            order=order,
            color="#2f6f5e" if perturbation_type == "canonical_watson_crick" else "#b2533e",
            size=5,
            jitter=0.12,
            alpha=0.85,
            ax=axis,
        )
        medians = panel.groupby("metric_label").value.median()
        for index, label in enumerate(order):
            if label in medians:
                axis.plot([index - 0.22, index + 0.22], [medians[label]] * 2, color="black", lw=2)
        axis.plot([2.65, 3.35], [0.5, 0.5], color="#666666", lw=1, ls="--")
        axis.set_title(title, fontsize=10)
        axis.set_xlabel("")
        axis.tick_params(axis="x", rotation=22)
        axis.set_ylim(-0.03, 1.03)
    axes[0].set_ylabel("Technical split agreement")
    axes[1].set_ylabel("")
    fig.suptitle("SaMBA spot-level technical repeatability", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGURES / "figure3_external_identifiability.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "figure3_external_identifiability.pdf", bbox_inches="tight")
    plt.close(fig)
    return source


def main() -> None:
    parser = argparse.ArgumentParser(description="Recover and evaluate Afek et al. SaMBA technical spots.")
    parser.add_argument(
        "--rebuild-spots",
        action="store_true",
        help="Reparse the official XLSX instead of reusing the processed Parquet cache.",
    )
    parser.add_argument(
        "--skip-secondary-splits",
        action="store_true",
        help="Skip the pre-registered 100-seed secondary split stability table.",
    )
    parser.add_argument(
        "--skip-permutations",
        action="store_true",
        help="Skip 2,000-permutation null tables.",
    )
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)

    if args.rebuild_spots or not SPOT_TABLE.exists():
        spots, _ = recover_spot_table()
    else:
        spots = pd.read_parquet(SPOT_TABLE)
    effects = split_effects(spots, SEED)
    effects.to_parquet(RESULTS / "external_samba_primary_split_effects.parquet", index=False)
    decomposition = primary_decomposition(effects)
    per_tf, _ = tf_level_summary(decomposition)
    render_figure3(per_tf)
    if not args.skip_secondary_splits:
        secondary_split_stability(spots)
    if not args.skip_permutations:
        permutation_nulls(effects, decomposition)
    print(
        f"SaMBA: {len(spots):,} spot rows, {len(effects):,} perturbations, "
        f"{decomposition.binding_site_id.nunique()} sites, {decomposition.tf_id.nunique()} TFs"
    )


if __name__ == "__main__":
    main()
