from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, gc_content, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
INTERIM_DIR = ensure_dir(ROOT / "data" / "interim")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")


BENCHMARK_COLUMNS = [
    "pair_id",
    "protein_name",
    "protein_sequence",
    "target_dna",
    "candidate_dna",
    "candidate_type",
    "mutation_count",
    "hamming_distance",
    "sequence_identity",
    "gc_content",
    "delta_gc",
    "source",
]


def load_validated_pairs() -> pd.DataFrame:
    validated_path = INTERIM_DIR / "dbp_target_pairs_validated.csv"
    if validated_path.exists():
        df = pd.read_csv(validated_path)
    else:
        raw_path = INTERIM_DIR / "dbp_target_pairs_raw.csv"
        if not raw_path.exists():
            raise FileNotFoundError("Missing raw or validated pair table")
        df = pd.read_csv(raw_path)
        df["keep_row"] = True
    df = df[df["keep_row"]].copy()
    df["protein_sequence"] = df["protein_sequence"].map(normalize_sequence)
    df["target_dna"] = df["target_dna"].map(normalize_sequence)
    df["dna_length"] = df["target_dna"].map(len)
    df["protein_length"] = df["protein_sequence"].map(len)
    df["gc_content"] = df["target_dna"].map(gc_content)
    return df.reset_index(drop=True)


def build_target_rows(pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in pairs.iterrows():
        rows.append(
            {
                "pair_id": row["pair_id"],
                "protein_name": row["protein_name"],
                "protein_sequence": row["protein_sequence"],
                "target_dna": row["target_dna"],
                "candidate_dna": row["target_dna"],
                "candidate_type": "target",
                "mutation_count": 0,
                "hamming_distance": 0,
                "sequence_identity": 1.0,
                "gc_content": float(row["gc_content"]),
                "delta_gc": 0.0,
                "source": row["source_name"],
            }
        )
    return pd.DataFrame(rows)


def build_benchmark() -> pd.DataFrame:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")
    target_rows = build_target_rows(pairs)
    parts = [target_rows]
    for filename in ["single_mutants.csv", "double_mutants.csv", "random_negatives.csv"]:
        path = PROCESSED_DIR / filename
        if path.exists():
            parts.append(pd.read_csv(path))
    benchmark = pd.concat(parts, ignore_index=True, sort=False)
    for col in BENCHMARK_COLUMNS:
        if col not in benchmark.columns:
            benchmark[col] = pd.NA
    benchmark = benchmark[BENCHMARK_COLUMNS + [c for c in benchmark.columns if c not in BENCHMARK_COLUMNS]]
    benchmark = benchmark.sort_values(["pair_id", "candidate_type", "mutation_count"], kind="stable").reset_index(drop=True)
    return benchmark


def main() -> None:
    pairs = load_validated_pairs()
    processed_pairs = pairs[
        [
            "pair_id",
            "pdb_id",
            "protein_entity_id",
            "protein_chain_ids",
            "protein_name",
            "protein_sequence",
            "protein_length",
            "dna_entity_id",
            "dna_chain_ids",
            "target_dna",
            "dna_length",
            "gc_content",
            "source_type",
            "source_name",
            "source_id",
            "source_url",
            "source_api_url",
            "source_title",
            "source_paper_title",
            "source_paper_doi",
            "source_paper_pmid",
            "source_paper_year",
            "retrieval_date",
            "experimental_or_designed",
            "has_specificity_ground_truth",
            "notes",
        ]
    ].copy()
    processed_pairs.to_csv(PROCESSED_DIR / "dbp_target_pairs.csv", index=False)

    # Build benchmark if mutant/control files exist.
    if (PROCESSED_DIR / "single_mutants.csv").exists() and (PROCESSED_DIR / "double_mutants.csv").exists():
        benchmark = build_benchmark()
        benchmark.to_csv(PROCESSED_DIR / "benchmark_v0_1.csv", index=False)
        summary = benchmark.groupby("candidate_type", dropna=False).size().reset_index(name="count")
        summary.to_csv(RESULTS_TABLES / "benchmark_summary.csv", index=False)
        print(summary.to_string(index=False))
        print(f"Benchmark written to {PROCESSED_DIR / 'benchmark_v0_1.csv'}")
    else:
        print(f"Processed pairs written to {PROCESSED_DIR / 'dbp_target_pairs.csv'}")


if __name__ == "__main__":
    main()

