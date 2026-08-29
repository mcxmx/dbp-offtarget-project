from __future__ import annotations

import random
import sys
import zlib
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import compute_sequence_metrics, ensure_dir, gc_content, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")
DNA_BASES = "ACGT"


def exact_gc_sequence(length: int, gc_count: int, rng: random.Random) -> str:
    if gc_count < 0 or gc_count > length:
        raise ValueError("gc_count out of range")
    gc_positions = set(rng.sample(range(length), gc_count)) if gc_count else set()
    seq = []
    for pos in range(length):
        if pos in gc_positions:
            seq.append(rng.choice("GC"))
        else:
            seq.append(rng.choice("AT"))
    return "".join(seq)


def best_gc_count(target_gc: float, length: int, tolerance: float) -> int:
    candidates = [k for k in range(length + 1) if abs(k / length - target_gc) <= tolerance]
    if not candidates:
        raise ValueError(f"No GC count within tolerance for length={length} and target_gc={target_gc}")
    return min(candidates, key=lambda k: abs(k / length - target_gc))


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(DNA_BASES) for _ in range(length))


def build_rows(pair: pd.Series, n_gc: int, n_random: int, seed: int, tolerance: float) -> list[dict]:
    target = normalize_sequence(pair["target_dna"])
    length = len(target)
    target_gc = gc_content(target)
    gc_count = best_gc_count(target_gc, length, tolerance)
    rows = []
    stable_offset = zlib.crc32(str(pair["pair_id"]).encode("utf-8")) % 10_000
    rng_gc = random.Random(seed + stable_offset)
    for i in range(n_gc):
        candidate = exact_gc_sequence(length, gc_count, rng_gc)
        metrics = compute_sequence_metrics(target, candidate, tuple(CONFIG["sequence_baseline_k_values"]))
        rows.append(
            {
                "pair_id": pair["pair_id"],
                "protein_name": pair["protein_name"],
                "protein_sequence": pair["protein_sequence"],
                "target_dna": target,
                "candidate_dna": candidate,
                "candidate_type": "gc_matched_random",
                "mutation_count": 0,
                "hamming_distance": metrics["hamming_distance"],
                "sequence_identity": metrics["sequence_identity"],
                "gc_content": metrics["gc_content"],
                "delta_gc": metrics["delta_gc"],
                "source": "generated_gc_matched",
                "target_gc_content": target_gc,
                "gc_tolerance": tolerance,
            }
        )
    rng_random = random.Random(seed + 17 + stable_offset)
    for i in range(n_random):
        candidate = random_dna(length, rng_random)
        metrics = compute_sequence_metrics(target, candidate, tuple(CONFIG["sequence_baseline_k_values"]))
        rows.append(
            {
                "pair_id": pair["pair_id"],
                "protein_name": pair["protein_name"],
                "protein_sequence": pair["protein_sequence"],
                "target_dna": target,
                "candidate_dna": candidate,
                "candidate_type": "random_dna",
                "mutation_count": 0,
                "hamming_distance": metrics["hamming_distance"],
                "sequence_identity": metrics["sequence_identity"],
                "gc_content": metrics["gc_content"],
                "delta_gc": metrics["delta_gc"],
                "source": "generated_random_dna",
                "target_gc_content": target_gc,
                "gc_tolerance": tolerance,
            }
        )
    return rows


def main() -> None:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")
    rows = []
    summary = []
    for _, pair in pairs.iterrows():
        pair_rows = build_rows(
            pair,
            int(CONFIG["gc_matched_per_target"]),
            int(CONFIG["random_per_target"]),
            int(CONFIG["seed"]),
            float(CONFIG["gc_tolerance"]),
        )
        rows.extend(pair_rows)
        summary.append(
            {
                "pair_id": pair["pair_id"],
                "dna_length": len(pair["target_dna"]),
                "target_gc_content": gc_content(pair["target_dna"]),
                "gc_matched_written": int(CONFIG["gc_matched_per_target"]),
                "random_written": int(CONFIG["random_per_target"]),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(PROCESSED_DIR / "random_negatives.csv", index=False)
    pd.DataFrame(summary).to_csv(RESULTS_TABLES / "random_controls_summary.csv", index=False)
    print(f"Random negatives: {len(out)}")


if __name__ == "__main__":
    main()
