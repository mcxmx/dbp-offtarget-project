from __future__ import annotations

from src.utils import is_valid_dna, normalize_sequence


DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    normalized = normalize_sequence(seq)
    if not is_valid_dna(normalized):
        raise ValueError(f"Invalid DNA sequence for reverse complement: {seq!r}")
    return normalized.translate(DNA_COMPLEMENT)[::-1]


def canonical_rc(seq: str) -> str:
    normalized = normalize_sequence(seq)
    rc = reverse_complement(normalized)
    return min(normalized, rc)


def rc_equivalent(seq1: str, seq2: str) -> bool:
    normalized_1 = normalize_sequence(seq1)
    normalized_2 = normalize_sequence(seq2)
    return canonical_rc(normalized_1) == canonical_rc(normalized_2)


def canonical_rc_set(sequences: list[str] | set[str] | tuple[str, ...]) -> set[str]:
    return {canonical_rc(sequence) for sequence in sequences}


def has_rc_split_leakage(train_sequences: list[str] | set[str], test_sequences: list[str] | set[str]) -> bool:
    return bool(canonical_rc_set(train_sequences) & canonical_rc_set(test_sequences))
