from __future__ import annotations

import math
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

import yaml

DNA_ALPHABET = "ACGT"
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")
DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_yaml(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def save_yaml(obj: dict, path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(obj, handle, sort_keys=False, allow_unicode=False)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def normalize_sequence(seq: str | None) -> str:
    if seq is None:
        return ""
    return re.sub(r"\s+", "", str(seq)).upper()


def gc_content(seq: str | None) -> float:
    seq = normalize_sequence(seq)
    if not seq:
        return math.nan
    gc = sum(base in {"G", "C"} for base in seq)
    return gc / len(seq)


def reverse_complement(seq: str | None) -> str:
    return normalize_sequence(seq).translate(DNA_COMPLEMENT)[::-1]


def is_valid_dna(seq: str | None) -> bool:
    seq = normalize_sequence(seq)
    return bool(seq) and all(base in DNA_ALPHABET for base in seq)


def is_valid_protein(seq: str | None) -> bool:
    seq = normalize_sequence(seq)
    return bool(seq) and all(residue in PROTEIN_ALPHABET for residue in seq)


def hamming_distance(a: str, b: str) -> int:
    a = normalize_sequence(a)
    b = normalize_sequence(b)
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(x != y for x, y in zip(a, b))


def edit_distance(a: str, b: str) -> int:
    a = normalize_sequence(a)
    b = normalize_sequence(b)
    if not a:
        return len(b)
    if not b:
        return len(a)
    if len(a) < len(b):
        a, b = b, a
    previous = list(range(len(b) + 1))
    for i, char_a in enumerate(a, start=1):
        current = [i]
        for j, char_b in enumerate(b, start=1):
            insert_cost = current[j - 1] + 1
            delete_cost = previous[j] + 1
            substitute_cost = previous[j - 1] + (char_a != char_b)
            current.append(min(insert_cost, delete_cost, substitute_cost))
        previous = current
    return previous[-1]


def sequence_identity(a: str, b: str) -> float:
    a = normalize_sequence(a)
    b = normalize_sequence(b)
    if not a and not b:
        return 1.0
    max_len = max(len(a), len(b))
    return max(0.0, 1.0 - edit_distance(a, b) / max_len)


def canonical_kmer_set(seq: str, k: int, rc_aware: bool = False) -> set[str]:
    seq = normalize_sequence(seq)
    if k <= 0 or len(seq) < k:
        return set()
    kmers = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if rc_aware:
            kmer = min(kmer, reverse_complement(kmer))
        kmers.add(kmer)
    return kmers


def kmer_jaccard(a: str, b: str, k: int, rc_aware: bool = False) -> float:
    set_a = canonical_kmer_set(a, k, rc_aware=rc_aware)
    set_b = canonical_kmer_set(b, k, rc_aware=rc_aware)
    if not set_a and not set_b:
        return 1.0
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


def gc_similarity(a: str, b: str) -> float:
    return max(0.0, 1.0 - abs(gc_content(a) - gc_content(b)))


def compute_sequence_metrics(target: str, candidate: str, k_values: tuple[int, ...] = (3, 4)) -> dict[str, float | int | str]:
    target = normalize_sequence(target)
    candidate = normalize_sequence(candidate)
    max_len = max(len(target), len(candidate), 1)
    min_len = min(len(target), len(candidate))
    hamming = hamming_distance(target, candidate) if len(target) == len(candidate) and target and candidate else math.nan
    edit = edit_distance(target, candidate)
    identity = sequence_identity(target, candidate)
    k3 = kmer_jaccard(target, candidate, 3, rc_aware=False)
    k4 = kmer_jaccard(target, candidate, 4, rc_aware=False)
    rck4 = kmer_jaccard(target, candidate, 4, rc_aware=True)
    gc_diff = gc_content(candidate) - gc_content(target)
    gc_sim = gc_similarity(target, candidate)
    proxy_score = (
        0.40 * identity
        + 0.15 * max(0.0, 1.0 - edit / max_len)
        + 0.20 * k3
        + 0.15 * k4
        + 0.10 * rck4
    )
    proxy_score = 0.9 * proxy_score + 0.1 * gc_sim
    return {
        "target_length": len(target),
        "candidate_length": len(candidate),
        "hamming_distance": hamming,
        "edit_distance": edit,
        "sequence_identity": identity,
        "gc_content": gc_content(candidate),
        "target_gc_content": gc_content(target),
        "delta_gc": gc_diff,
        "kmer3_jaccard": k3,
        "kmer4_jaccard": k4,
        "rc_kmer4_jaccard": rck4,
        "gc_similarity": gc_sim,
        "proxy_score": proxy_score,
    }


def count_kmers(seq: str, k: int, rc_aware: bool = False) -> Counter:
    seq = normalize_sequence(seq)
    counts: Counter = Counter()
    if k <= 0 or len(seq) < k:
        return counts
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if rc_aware:
            kmer = min(kmer, reverse_complement(kmer))
        counts[kmer] += 1
    return counts


def clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def canonical_dna(seq: str | None) -> str:
    return normalize_sequence(seq)


def reverse_complement_canonical(seq: str | None) -> str:
    seq = normalize_sequence(seq)
    rc = reverse_complement(seq)
    return min(seq, rc)
