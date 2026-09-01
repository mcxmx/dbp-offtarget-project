from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import compute_sequence_metrics, ensure_dir, load_yaml, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
INTERIM_DIR = ensure_dir(ROOT / "data" / "interim" / "gse237017")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed" / "v0_3")
METADATA_DIR = ensure_dir(ROOT / "metadata" / "v0_3")
TABLES_DIR = ensure_dir(ROOT / "results" / "v0_3" / "tables")


def target_7mers(target: str) -> list[str]:
    if not isinstance(target, str) or len(target) < 7:
        return []
    return [target[i : i + 7] for i in range(len(target) - 6)]


def spearman(x: pd.Series, y: pd.Series) -> float:
    return x.rank(method="average").corr(y.rank(method="average"), method="pearson")


def max_similarity_to_targets(sequence: str, refs: list[str]) -> dict[str, float | str]:
    if not refs:
        return {
            "best_matching_target_7mer": "",
            "hamming_similarity_to_target_7mer": np.nan,
            "edit_similarity_to_target_7mer": np.nan,
            "kmer3_jaccard_to_target_7mer": np.nan,
            "kmer4_jaccard_to_target_7mer": np.nan,
        }
    best = None
    for ref in refs:
        metrics = compute_sequence_metrics(ref, sequence)
        length = max(len(ref), 1)
        row = {
            "best_matching_target_7mer": ref,
            "hamming_similarity_to_target_7mer": 1.0 - metrics["hamming_distance"] / length,
            "edit_similarity_to_target_7mer": 1.0 - metrics["edit_distance"] / length,
            "kmer3_jaccard_to_target_7mer": metrics["kmer3_jaccard"],
            "kmer4_jaccard_to_target_7mer": metrics["kmer4_jaccard"],
        }
        if best is None or row["hamming_similarity_to_target_7mer"] > best["hamming_similarity_to_target_7mer"]:
            best = row
    return best or {}


def build_consensus(long_df: pd.DataFrame, sequences: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        long_df.groupby(["protein_id", "dna_7mer"], as_index=False)
        .agg(
            e_score_mean=("e_score", "mean"),
            e_score_median=("e_score", "median"),
            e_score_std=("e_score", "std"),
            n_replicates=("gsm_id", "nunique"),
            median_intensity_mean=("median_intensity", "mean"),
            median_intensity_median=("median_intensity", "median"),
            median_intensity_std=("median_intensity", "std"),
            z_score_mean=("z_score", "mean"),
            z_score_median=("z_score", "median"),
            z_score_std=("z_score", "std"),
            source_gsms=("gsm_id", lambda values: ";".join(sorted(set(values)))),
            source_files=("source_file", lambda values: ";".join(sorted(set(values)))),
            source_file_sha256=("source_file_sha256", lambda values: ";".join(sorted(set(values)))),
            protein_concentration=("protein_concentration", lambda values: ";".join(sorted(set(values)))),
        )
    )
    grouped["experimental_score_primary"] = grouped["e_score_mean"]
    grouped["experimental_score_type"] = "PBM E-score"
    grouped["experimental_score_raw"] = grouped["experimental_score_primary"]
    grouped["source_gse"] = CONFIG["benchmark_v0_3"]["gse_id"]
    grouped = grouped.merge(sequences[["protein_id", "protein_sequence", "sequence_length"]], on="protein_id", how="left")
    grouped = grouped.merge(
        targets[["protein_id", "intended_target_dna", "target_length", "target_id", "target_context", "target_duplex"]],
        on="protein_id",
        how="left",
    )
    normalized = []
    for _, group in grouped.groupby("protein_id"):
        out = group.copy()
        values = out["experimental_score_primary"]
        out["experimental_rank"] = values.rank(method="average", ascending=False)
        out["experimental_percentile"] = values.rank(method="average", pct=True)
        std = values.std(ddof=0)
        out["experimental_score_within_protein_z"] = (values - values.mean()) / std if std and not np.isnan(std) else np.nan
        normalized.append(out)
    return pd.concat(normalized, ignore_index=True).sort_values(["protein_id", "dna_7mer"]).reset_index(drop=True)


def protein_specificity_summary(benchmark: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein_id, group in benchmark.groupby("protein_id"):
        score = group["experimental_score_primary"]
        n = len(group)
        top_1_count = max(1, math.ceil(0.01 * n))
        top_5_count = max(1, math.ceil(0.05 * n))
        sorted_scores = score.sort_values(ascending=False)
        shifted = score - score.min()
        total_shifted = shifted.sum()
        rows.append(
            {
                "protein_id": protein_id,
                "n_7mers": n,
                "n_replicates_min": int(group["n_replicates"].min()),
                "n_replicates_max": int(group["n_replicates"].max()),
                "max_e_score": float(score.max()),
                "median_e_score": float(score.median()),
                "q95_e_score": float(score.quantile(0.95)),
                "q99_e_score": float(score.quantile(0.99)),
                "n_top_1_percent_7mers": top_1_count,
                "n_top_5_percent_7mers": top_5_count,
                "n_high_score_7mers": pd.NA,
                "high_score_threshold_source": "not_applied; no PBM E-score threshold was assumed for v0.3",
                "top1_percent_score_mass_fraction_derived": float(shifted.loc[sorted_scores.head(top_1_count).index].sum() / total_shifted) if total_shifted else np.nan,
                "specificity_entropy_or_concentration_metric": "top1_percent_score_mass_fraction_derived",
                "notes": "E-score distribution summarized per protein; absolute values are not treated as cross-protein binding affinity.",
            }
        )
    return pd.DataFrame(rows).sort_values("protein_id")


def target_rank_summary(benchmark: pd.DataFrame, targets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    detail_rows = []
    for _, target_row in targets.sort_values("protein_id").iterrows():
        protein_id = target_row["protein_id"]
        target = target_row["intended_target_dna"]
        kmers = target_7mers(target)
        group = benchmark[benchmark["protein_id"] == protein_id].set_index("dna_7mer")
        target_scores = []
        for kmer in kmers:
            if kmer in group.index:
                row = group.loc[kmer]
                target_scores.append(float(row["experimental_score_primary"]))
                detail_rows.append(
                    {
                        "protein_id": protein_id,
                        "intended_target_dna": target,
                        "target_7mer": kmer,
                        "e_score": float(row["experimental_score_primary"]),
                        "experimental_rank": float(row["experimental_rank"]),
                        "experimental_percentile": float(row["experimental_percentile"]),
                        "notes": "Target-derived 7-mer score; not a full-target binding score.",
                    }
                )
        if target_scores:
            detail = pd.DataFrame([row for row in detail_rows if row["protein_id"] == protein_id])
            best_idx = detail["e_score"].idxmax()
            best = detail.loc[best_idx]
            summary_rows.append(
                {
                    "protein_id": protein_id,
                    "target_length": int(target_row["target_length"]),
                    "n_target_7mers": len(kmers),
                    "best_target_7mer": best["target_7mer"],
                    "best_target_7mer_escore": float(best["e_score"]),
                    "best_target_percentile": float(best["experimental_percentile"]),
                    "mean_target_7mer_escore": float(np.mean(target_scores)),
                    "median_target_7mer_escore": float(np.median(target_scores)),
                    "notes": "Intended target length exceeds 7 bp for all designs; summary uses overlapping target-derived 7-mers and is not a full-target affinity.",
                }
            )
        else:
            summary_rows.append(
                {
                    "protein_id": protein_id,
                    "target_length": int(target_row["target_length"]),
                    "n_target_7mers": len(kmers),
                    "best_target_7mer": "",
                    "best_target_7mer_escore": np.nan,
                    "best_target_percentile": np.nan,
                    "mean_target_7mer_escore": np.nan,
                    "median_target_7mer_escore": np.nan,
                    "notes": "No target-derived 7-mer found in benchmark table.",
                }
            )
    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


def sequence_baseline(benchmark: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scored_rows = []
    for protein_id, group in benchmark.groupby("protein_id"):
        target = group["intended_target_dna"].iloc[0]
        refs = target_7mers(target)
        for _, row in group.iterrows():
            scored_rows.append({**row[["protein_id", "dna_7mer", "experimental_score_primary"]].to_dict(), **max_similarity_to_targets(row["dna_7mer"], refs)})
    scored = pd.DataFrame(scored_rows)
    metrics = [
        "hamming_similarity_to_target_7mer",
        "edit_similarity_to_target_7mer",
        "kmer3_jaccard_to_target_7mer",
        "kmer4_jaccard_to_target_7mer",
    ]
    rows = []
    for (protein_id, group) in scored.groupby("protein_id"):
        for metric in metrics:
            rows.append(
                {
                    "protein_id": protein_id,
                    "metric": metric,
                    "spearman": spearman(group[metric], group["experimental_score_primary"]),
                    "n_sequences": int(len(group)),
                    "reference_definition": "max similarity to any overlapping 7-mer from the intended target DNA",
                    "notes": "Per-protein sequence-only baseline; not protein-conditioned prediction.",
                }
            )
    return pd.DataFrame(rows), scored


def main() -> None:
    long_df = pd.read_parquet(INTERIM_DIR / "upbm_7mers_long.parquet")
    sequences = pd.read_csv(METADATA_DIR / "designed_dbp_sequences.csv")
    targets = pd.read_csv(METADATA_DIR / "designed_dbp_targets.csv")
    benchmark = build_consensus(long_df, sequences, targets)
    benchmark.to_parquet(PROCESSED_DIR / "designed_dbp_upbm_v0_3.parquet", index=False)
    benchmark.head(200).to_csv(PROCESSED_DIR / "designed_dbp_upbm_v0_3_preview.csv", index=False)
    protein_specificity_summary(benchmark).to_csv(TABLES_DIR / "protein_specificity_summary.csv", index=False)
    target_summary, target_detail = target_rank_summary(benchmark, targets)
    target_summary.to_csv(TABLES_DIR / "target_rank_summary.csv", index=False)
    target_detail.to_csv(TABLES_DIR / "target_7mer_scores.csv", index=False)
    baseline, baseline_scored = sequence_baseline(benchmark)
    baseline.to_csv(TABLES_DIR / "designed_dbp_sequence_baseline.csv", index=False)
    baseline_scored.to_parquet(PROCESSED_DIR / "designed_dbp_sequence_baseline_scored_v0_3.parquet", index=False)
    print(f"benchmark rows: {len(benchmark)}")
    print(f"proteins: {benchmark['protein_id'].nunique()}")
    print(benchmark.groupby("protein_id")["dna_7mer"].nunique().to_string())
    print(target_summary.to_string(index=False))
    print(baseline.to_string(index=False))


if __name__ == "__main__":
    main()
