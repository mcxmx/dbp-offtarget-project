from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.models.base_binding_model import SequenceProxyBaseline
from src.utils import compute_sequence_metrics, ensure_dir, load_yaml, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")


def score_benchmark(df: pd.DataFrame) -> pd.DataFrame:
    scored_rows = []
    proxy = SequenceProxyBaseline(tuple(CONFIG["sequence_baseline_k_values"]), CONFIG["proxy_weights"])
    for _, row in df.iterrows():
        metrics = compute_sequence_metrics(row["target_dna"], row["candidate_dna"], tuple(CONFIG["sequence_baseline_k_values"]))
        scored_rows.append(
            {
                **row.to_dict(),
                **metrics,
                "sequence_proxy_score": metrics["proxy_score"],
                "sequence_proxy_model_score": proxy.score(row["protein_sequence"], row["candidate_dna"]),
                "sequence_proxy_model_is_protein_conditioned": proxy.is_protein_conditioned,
                "sequence_proxy_model_label": proxy.score_label,
            }
        )
    return pd.DataFrame(scored_rows)


def main() -> None:
    benchmark_path = PROCESSED_DIR / "benchmark_v0_1.csv"
    if not benchmark_path.exists():
        raise FileNotFoundError(f"Missing benchmark table: {benchmark_path}")
    df = pd.read_csv(benchmark_path)
    scored = score_benchmark(df)
    scored.to_csv(PROCESSED_DIR / "benchmark_v0_1_scored.csv", index=False)

    summary = (
        scored.groupby("candidate_type")
        .agg(
            n=("pair_id", "size"),
            mean_proxy_score=("sequence_proxy_score", "mean"),
            std_proxy_score=("sequence_proxy_score", "std"),
            mean_gc_content=("gc_content", "mean"),
            mean_identity=("sequence_identity", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(RESULTS_TABLES / "sequence_baseline_summary.csv", index=False)
    scored[[
        "pair_id",
        "candidate_type",
        "sequence_proxy_score",
        "sequence_proxy_model_score",
        "sequence_proxy_model_is_protein_conditioned",
        "sequence_proxy_model_label",
        "hamming_distance",
        "edit_distance",
        "sequence_identity",
        "kmer3_jaccard",
        "kmer4_jaccard",
        "rc_kmer4_jaccard",
        "gc_similarity",
    ]].to_csv(RESULTS_TABLES / "sequence_baseline_metrics.csv", index=False)

    print(summary.to_string(index=False))
    print(f"Scored benchmark written to {PROCESSED_DIR / 'benchmark_v0_1_scored.csv'}")


if __name__ == "__main__":
    main()
