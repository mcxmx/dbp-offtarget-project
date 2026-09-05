from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import numpy as np
import pandas as pd

from src.v0_5_training import sample_rank_pairs


EASY_MIN_RANK_DIFFERENCE = 0.50
MEDIUM_MIN_RANK_DIFFERENCE = 0.20
N_RANK_DECILES = 10


def rank_deciles(scores: np.ndarray) -> np.ndarray:
    ranks = pd.Series(scores).rank(method="average", pct=True).to_numpy()
    return np.minimum((ranks * N_RANK_DECILES).astype(int), N_RANK_DECILES - 1)


def _difficulty(rank_difference: float) -> str | None:
    if rank_difference >= EASY_MIN_RANK_DIFFERENCE:
        return "easy"
    if rank_difference >= MEDIUM_MIN_RANK_DIFFERENCE:
        return "medium"
    if rank_difference > 0.0:
        return "hard"
    return None


def _target_counts(pair_count: int) -> dict[str, int]:
    easy = int(round(pair_count * 0.40))
    medium = int(round(pair_count * 0.35))
    hard = pair_count - easy - medium
    return {"easy": easy, "medium": medium, "hard": hard}


def _stable_seed(protein: str, seed: int, protocol: str) -> int:
    text = f"{protein}|{protocol}"
    return seed + sum((index + 1) * ord(char) for index, char in enumerate(text)) % 1_000_003


def _sample_dense_for_protein(
    group: pd.DataFrame,
    *,
    protein: str,
    pair_count: int,
    seed: int,
    protocol: str,
    tie_tolerance: float,
    min_coverage: float,
) -> pd.DataFrame:
    group = group.reset_index(drop=True)
    scores = group["experimental_score"].to_numpy(dtype=float)
    ranks = pd.Series(scores).rank(method="average", pct=True).to_numpy()
    deciles = rank_deciles(scores)
    rng = np.random.default_rng(_stable_seed(protein, seed, protocol))
    counts = _target_counts(pair_count)
    selected_counts = {"easy": 0, "medium": 0, "hard": 0}
    selected: list[dict[str, Any]] = []
    selected_keys: set[tuple[int, int]] = set()
    covered: set[int] = set()

    def try_add(label: str, left: int, right: int) -> bool:
        if left == right or abs(scores[left] - scores[right]) <= tie_tolerance:
            return False
        if _difficulty(abs(ranks[left] - ranks[right])) != label:
            return False
        key = tuple(sorted((left, right)))
        if key in selected_keys:
            return False
        selected_keys.add(key)
        selected.append(
            {
                "protein_id": protein,
                "left_index": left,
                "right_index": right,
                "left_score": scores[left],
                "right_score": scores[right],
                "difficulty": label,
            }
        )
        selected_counts[label] += 1
        covered.update(key)
        return True

    rank_order = np.argsort(ranks, kind="mergesort")
    offset_fractions = {
        "easy": (0.51, 0.60, 0.70, 0.80, 0.90),
        "medium": (0.21, 0.30, 0.40, 0.49),
        "hard": tuple(index / max(len(group), 1) for index in range(1, min(32, len(group)))),
    }

    def candidate_pairs(label: str) -> list[tuple[int, int]]:
        pairs: list[tuple[int, int]] = []
        seen: set[tuple[int, int]] = set()
        permutation = rng.permutation(len(group))
        for fraction in offset_fractions[label]:
            offset = max(1, int(np.ceil(len(group) * fraction)))
            for position in permutation:
                for direction in (1, -1):
                    partner_position = int(position) + direction * offset
                    if partner_position < 0 or partner_position >= len(group):
                        continue
                    left = int(rank_order[int(position)])
                    right = int(rank_order[partner_position])
                    if _difficulty(abs(ranks[left] - ranks[right])) != label:
                        continue
                    key = tuple(sorted((left, right)))
                    if key not in seen:
                        seen.add(key)
                        pairs.append((left, right))
        return pairs

    # Construct rank-ordered candidate pools rather than repeatedly rejecting
    # random pairs. The first pass prioritizes unseen endpoints; the second
    # pass fills the exact difficulty quota.
    for label, count in counts.items():
        candidates = candidate_pairs(label)
        target_endpoints = int(np.ceil(len(group) * min_coverage))
        for left, right in candidates:
            if selected_counts[label] >= count:
                break
            if left not in covered or right not in covered:
                try_add(label, left, right)
        for left, right in candidates:
            if selected_counts[label] >= count:
                break
            try_add(label, left, right)
        if selected_counts[label] < count:
            raise RuntimeError(
                f"Unable to sample {count} unique {label} pairs for {protein}; "
                f"selected={selected_counts[label]}, "
                f"candidate_pool={len(candidates)}"
            )
        if len(covered) < target_endpoints:
            # Continue with other strata before declaring the registered
            # protocol target unreachable; coverage is audited after sampling.
            continue

    return pd.DataFrame(selected)


def sample_dense_pairs(
    benchmark: pd.DataFrame,
    proteins: Iterable[str],
    *,
    pairs_per_protein: int,
    seed: int,
    protocol: str,
    tie_tolerance: float = 1e-10,
) -> pd.DataFrame:
    if protocol not in {"D4096", "D16384"}:
        raise ValueError(f"Dense sampler requires D4096 or D16384, got {protocol}")
    if pairs_per_protein <= 0:
        raise ValueError("pairs_per_protein must be positive")
    min_coverage = 0.60 if protocol == "D4096" else 0.90
    parts = []
    for protein in sorted(set(proteins)):
        group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
        if group.empty:
            raise ValueError(f"No benchmark rows for {protein}")
        parts.append(
            _sample_dense_for_protein(
                group,
                protein=protein,
                pair_count=pairs_per_protein,
                seed=seed,
                protocol=protocol,
                tie_tolerance=tie_tolerance,
                min_coverage=min_coverage,
            )
        )
    return pd.concat(parts, ignore_index=True)


def audit_pair_sampling(
    benchmark: pd.DataFrame,
    proteins: Iterable[str],
    *,
    seed: int = 17,
    tie_tolerance: float = 1e-10,
) -> pd.DataFrame:
    rows = []
    for protocol, pair_count in [("S512", 512), ("D4096", 4096), ("D16384", 16384)]:
        if protocol == "S512":
            pairs = sample_rank_pairs(
                benchmark,
                list(proteins),
                pairs_per_protein=pair_count,
                seed=seed,
                tie_tolerance=tie_tolerance,
            )
        else:
            pairs = sample_dense_pairs(
                benchmark,
                proteins,
                pairs_per_protein=pair_count,
                seed=seed,
                protocol=protocol,
                tie_tolerance=tie_tolerance,
            )
        for protein in sorted(set(proteins)):
            group = benchmark.loc[benchmark["protein_id"].eq(protein)].reset_index(drop=True)
            protein_pairs = pairs.loc[pairs["protein_id"].eq(protein)].copy()
            endpoints = set(protein_pairs["left_index"]) | set(protein_pairs["right_index"])
            scores = group["experimental_score"].to_numpy(dtype=float)
            deciles = rank_deciles(scores)
            covered_deciles = sorted({int(deciles[index]) for index in endpoints})
            pair_keys = protein_pairs.apply(
                lambda row: tuple(sorted((int(row["left_index"]), int(row["right_index"])))),
                axis=1,
            )
            rows.append(
                {
                    "protein": protein,
                    "protocol": protocol,
                    "pair_count": len(protein_pairs),
                    "unique_candidates": len(endpoints),
                    "candidate_coverage": len(endpoints) / len(group),
                    "rank_deciles_covered": "|".join(map(str, covered_deciles)),
                    "n_rank_deciles_covered": len(covered_deciles),
                    "easy_pairs": int((protein_pairs["difficulty"] == "easy").sum()),
                    "medium_pairs": int((protein_pairs["difficulty"] == "medium").sum()),
                    "hard_pairs": int((protein_pairs["difficulty"] == "hard").sum()),
                    "duplicate_pair_fraction": 1.0 - pair_keys.nunique() / len(protein_pairs),
                }
            )
    return pd.DataFrame(rows)
