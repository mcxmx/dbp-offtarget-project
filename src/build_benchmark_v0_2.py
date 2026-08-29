from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.generate_gc_matched import build_rows as random_control_rows
from src.generate_mutants import double_mutant_rows, maybe_sample, single_mutant_rows
from src.utils import compute_sequence_metrics, ensure_dir, gc_content, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")

PAIR_COLUMNS = [
    "source_type",
    "source_name",
    "source_id",
    "source_url",
    "source_title",
    "source_paper_title",
    "source_paper_doi",
    "source_paper_pmid",
    "retrieval_date",
    "uniprot_id",
    "dna_role",
    "binding_mechanism",
    "sequence_specificity_class",
    "recommended_use",
    "has_structural_cognate",
    "has_direct_dna_binding_evidence",
    "has_sequence_specificity_evidence",
    "has_quantitative_specificity_ground_truth",
    "curation_confidence",
    "curation_reason",
]


def pair_metadata(pairs: pd.DataFrame) -> pd.DataFrame:
    keep = ["pair_id"] + [col for col in PAIR_COLUMNS if col in pairs.columns]
    return pairs[keep].copy()


def add_pair_metadata(df: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    return df.merge(pair_metadata(pairs), on="pair_id", how="left")


def target_rows(pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, pair in pairs.iterrows():
        target = normalize_sequence(pair["target_dna"])
        rows.append(
            {
                "pair_id": pair["pair_id"],
                "protein_name": pair["protein_name"],
                "protein_sequence": pair["protein_sequence"],
                "target_dna": target,
                "candidate_dna": target,
                "candidate_type": "target",
                "mutation_type": "none",
                "mutation_count": 0,
                "mutation_positions": "",
                "original_bases": "",
                "mutated_bases": "",
                "hamming_distance": 0,
                "edit_distance": 0,
                "sequence_identity": 1.0,
                "gc_content": gc_content(target),
                "delta_gc": 0.0,
                "source": "curated_structural_cognate",
            }
        )
    return pd.DataFrame(rows)


def score_benchmark(df: pd.DataFrame) -> pd.DataFrame:
    metric_rows = []
    for _, row in df.iterrows():
        metrics = compute_sequence_metrics(
            row["target_dna"],
            row["candidate_dna"],
            tuple(CONFIG["sequence_baseline_k_values"]),
        )
        metric_rows.append(metrics)
    metric_df = pd.DataFrame(metric_rows)
    out = df.copy()
    for col in metric_df.columns:
        out[col] = metric_df[col]
    out = out.rename(columns={"proxy_score": "sequence_proxy_score"})
    return out


def duplicate_qc(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (pair_id, candidate_type), group in df.groupby(["pair_id", "candidate_type"]):
        rows.append(
            {
                "pair_id": pair_id,
                "candidate_type": candidate_type,
                "rows": len(group),
                "unique_candidate_dna": group["candidate_dna"].nunique(),
                "duplicate_candidate_dna_rows": len(group) - group["candidate_dna"].nunique(),
                "duplicates_removed": False,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs_v0_2.csv")
    if pairs.empty:
        raise ValueError("No curated v0.2 pairs found")

    forbidden = set(CONFIG.get("benchmark_v0_2", {}).get("exclude_sequence_specificity_class", []))
    observed_forbidden = sorted(set(pairs["sequence_specificity_class"]) & forbidden)
    if observed_forbidden:
        raise AssertionError(f"Forbidden specificity classes in v0.2 core pairs: {observed_forbidden}")

    single_all = []
    double_all = []
    mutant_summary = []
    max_double = int(CONFIG["max_double_mutants_per_target"])
    for _, pair in pairs.iterrows():
        single_rows = single_mutant_rows(pair)
        double_rows = double_mutant_rows(pair)
        target_length = len(normalize_sequence(pair["target_dna"]))
        theoretical_single = 3 * target_length
        theoretical_double = 9 * math.comb(target_length, 2)
        if len(single_rows) != theoretical_single:
            raise AssertionError(f"{pair['pair_id']}: expected {theoretical_single} single mutants, got {len(single_rows)}")
        if len(double_rows) != theoretical_double:
            raise AssertionError(f"{pair['pair_id']}: expected {theoretical_double} double mutants, got {len(double_rows)}")
        total_double = len(double_rows)
        if len(double_rows) > max_double:
            double_rows, total_double = maybe_sample(double_rows, max_double, CONFIG["seed"])
        single_all.extend(single_rows)
        double_all.extend(double_rows)
        mutant_summary.append(
            {
                "pair_id": pair["pair_id"],
                "dna_length": target_length,
                "single_mutants_theoretical": theoretical_single,
                "single_mutants_written": len(single_rows),
                "double_mutants_theoretical": theoretical_double,
                "double_mutants_before_sampling": total_double,
                "double_mutants_written": len(double_rows),
                "double_sampling_applied": len(double_rows) < theoretical_double,
            }
        )

    single_df = add_pair_metadata(pd.DataFrame(single_all), pairs)
    double_df = add_pair_metadata(pd.DataFrame(double_all), pairs)

    random_rows = []
    random_summary = []
    for _, pair in pairs.iterrows():
        rows = random_control_rows(
            pair,
            int(CONFIG["gc_matched_per_target"]),
            int(CONFIG["random_per_target"]),
            int(CONFIG["seed"]),
            float(CONFIG["gc_tolerance"]),
        )
        random_rows.extend(rows)
        random_summary.append(
            {
                "pair_id": pair["pair_id"],
                "dna_length": len(normalize_sequence(pair["target_dna"])),
                "target_gc_content": gc_content(pair["target_dna"]),
                "gc_matched_written": int(CONFIG["gc_matched_per_target"]),
                "random_written": int(CONFIG["random_per_target"]),
            }
        )
    random_df = add_pair_metadata(pd.DataFrame(random_rows), pairs)
    target_df = add_pair_metadata(target_rows(pairs), pairs)

    benchmark = pd.concat([target_df, single_df, double_df, random_df], ignore_index=True, sort=False)
    scored = score_benchmark(benchmark)
    summary = (
        scored.groupby("candidate_type")
        .size()
        .reset_index(name="count")
        .sort_values("candidate_type")
    )
    summary.loc[len(summary)] = {"candidate_type": "total", "count": int(len(scored))}

    single_df.to_csv(PROCESSED_DIR / "single_mutants_v0_2.csv", index=False)
    double_df.to_csv(PROCESSED_DIR / "double_mutants_v0_2.csv", index=False)
    random_df.to_csv(PROCESSED_DIR / "random_negatives_v0_2.csv", index=False)
    benchmark.to_csv(PROCESSED_DIR / "benchmark_v0_2.csv", index=False)
    scored.to_csv(PROCESSED_DIR / "benchmark_v0_2_scored.csv", index=False)
    pd.DataFrame(mutant_summary).to_csv(RESULTS_TABLES / "mutant_counts_v0_2.csv", index=False)
    pd.DataFrame(random_summary).to_csv(RESULTS_TABLES / "random_controls_summary_v0_2.csv", index=False)
    duplicate_qc(benchmark).to_csv(RESULTS_TABLES / "benchmark_v0_2_generation_qc.csv", index=False)
    summary.to_csv(RESULTS_TABLES / "benchmark_summary_v0_2.csv", index=False)

    print(f"v0.2 pairs: {len(pairs)}")
    print(f"single mutants: {len(single_df)}")
    print(f"double mutants: {len(double_df)}")
    print(f"random negatives: {len(random_df)}")
    print(f"benchmark rows: {len(benchmark)}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

