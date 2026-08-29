from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, gc_content, load_yaml, normalize_sequence, project_root, compute_sequence_metrics


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")
DNA_BASES = "ACGT"


def single_mutant_rows(pair: pd.Series) -> list[dict]:
    target = normalize_sequence(pair["target_dna"])
    rows = []
    for pos, original in enumerate(target):
        for alt in DNA_BASES:
            if alt == original:
                continue
            candidate = target[:pos] + alt + target[pos + 1 :]
            metrics = compute_sequence_metrics(target, candidate, tuple(CONFIG["sequence_baseline_k_values"]))
            rows.append(
                {
                    "pair_id": pair["pair_id"],
                    "protein_name": pair["protein_name"],
                    "protein_sequence": pair["protein_sequence"],
                    "target_dna": target,
                    "candidate_dna": candidate,
                    "candidate_type": "single_mut",
                    "mutation_type": "single_mut",
                    "mutation_count": 1,
                    "mutation_positions": str(pos + 1),
                    "original_bases": original,
                    "mutated_bases": alt,
                    "hamming_distance": 1,
                    "edit_distance": metrics["edit_distance"],
                    "sequence_identity": metrics["sequence_identity"],
                    "gc_content": metrics["gc_content"],
                    "delta_gc": metrics["delta_gc"],
                    "source": "generated_single_mutant",
                }
            )
    return rows


def double_mutant_rows(pair: pd.Series) -> list[dict]:
    target = normalize_sequence(pair["target_dna"])
    rows = []
    positions = list(range(len(target)))
    for i, j in itertools.combinations(positions, 2):
        for alt_i in DNA_BASES:
            if alt_i == target[i]:
                continue
            for alt_j in DNA_BASES:
                if alt_j == target[j]:
                    continue
                candidate = list(target)
                candidate[i] = alt_i
                candidate[j] = alt_j
                candidate = "".join(candidate)
                metrics = compute_sequence_metrics(target, candidate, tuple(CONFIG["sequence_baseline_k_values"]))
                rows.append(
                    {
                        "pair_id": pair["pair_id"],
                        "protein_name": pair["protein_name"],
                        "protein_sequence": pair["protein_sequence"],
                        "target_dna": target,
                        "candidate_dna": candidate,
                        "candidate_type": "double_mut",
                        "mutation_type": "double_mut",
                        "mutation_count": 2,
                        "mutation_positions": f"{i + 1};{j + 1}",
                        "original_bases": f"{target[i]};{target[j]}",
                        "mutated_bases": f"{alt_i};{alt_j}",
                        "hamming_distance": 2,
                        "edit_distance": metrics["edit_distance"],
                        "sequence_identity": metrics["sequence_identity"],
                        "gc_content": metrics["gc_content"],
                        "delta_gc": metrics["delta_gc"],
                        "source": "generated_double_mutant",
                    }
                )
    return rows


def maybe_sample(rows: list[dict], max_rows: int, seed: int) -> tuple[list[dict], int]:
    if len(rows) <= max_rows:
        return rows, len(rows)
    sampled = pd.DataFrame(rows).sample(n=max_rows, random_state=seed).to_dict(orient="records")
    return sampled, len(rows)


def main() -> None:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")
    single_all = []
    double_all = []
    summary_rows = []
    for _, pair in pairs.iterrows():
        single_rows = single_mutant_rows(pair)
        double_rows = double_mutant_rows(pair)
        theoretical_single = 3 * len(pair["target_dna"])
        theoretical_double = 9 * math.comb(len(pair["target_dna"]), 2)
        if len(single_rows) != theoretical_single:
            raise AssertionError(f"{pair['pair_id']}: expected {theoretical_single} single mutants, got {len(single_rows)}")
        if len(double_rows) != theoretical_double:
            raise AssertionError(f"{pair['pair_id']}: expected {theoretical_double} double mutants, got {len(double_rows)}")
        max_double = int(CONFIG["max_double_mutants_per_target"])
        if len(double_rows) > max_double:
            double_rows, total_double = maybe_sample(double_rows, max_double, CONFIG["seed"])
        else:
            total_double = len(double_rows)
        single_all.extend(single_rows)
        double_all.extend(double_rows)
        summary_rows.append(
            {
                "pair_id": pair["pair_id"],
                "dna_length": len(pair["target_dna"]),
                "single_mutants_theoretical": theoretical_single,
                "double_mutants_theoretical": theoretical_double,
                "double_mutants_written": len(double_rows),
            }
        )

    single_df = pd.DataFrame(single_all)
    double_df = pd.DataFrame(double_all)
    single_df.to_csv(PROCESSED_DIR / "single_mutants.csv", index=False)
    double_df.to_csv(PROCESSED_DIR / "double_mutants.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(RESULTS_TABLES / "mutant_counts.csv", index=False)

    print(f"Single mutants: {len(single_df)}")
    print(f"Double mutants: {len(double_df)}")


if __name__ == "__main__":
    main()

