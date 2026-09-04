"""Project-side parser and fixed scorer for official DeepPBS predictions.

DeepPBS is an external, structure-aware baseline.  This module does not
reimplement its model; it validates the official NPZ output and maps its
position-wise base probabilities to the project's RC-class 7-mer units using
the frozen v0.4.2 protocol.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from src.sequence_equivalence import canonical_rc, reverse_complement


DEEPPBS_BASE_ORDER = ("A", "C", "G", "T")
BASE_TO_INDEX = {base: index for index, base in enumerate(DEEPPBS_BASE_ORDER)}
PRIMARY_EPSILON = 1e-9


def _one_hot_to_sequence(seq_array: np.ndarray) -> str:
    array = np.asarray(seq_array)
    if array.ndim == 1:
        if not np.all(np.isin(array, [0, 1])):
            raise ValueError("Seq vector must contain binary one-hot values")
        if array.size % 4 != 0:
            raise ValueError(f"Seq vector length is not divisible by 4: {array.shape}")
        array = array.reshape((-1, 4))
    if array.ndim != 2 or array.shape[1] != 4:
        raise ValueError(f"Seq must have shape (L, 4), got {array.shape}")
    if not np.allclose(array.sum(axis=1), 1.0):
        raise ValueError("Seq rows are not one-hot")
    indices = np.argmax(array, axis=1)
    return "".join(DEEPPBS_BASE_ORDER[int(index)] for index in indices)


def parse_prediction_npz(path: str | Path) -> dict[str, object]:
    """Load and validate a DeepPBS prediction NPZ.

    The upstream ``predict.py`` writes ``P`` after softmax, ensemble averaging,
    and strand averaging.  ``Seq`` is the hard input sequence tensor.  The
    caller still records the upstream source files in run provenance.
    """

    path = Path(path)
    with np.load(path, allow_pickle=False) as data:
        keys = set(data.files)
        required = {"P", "Seq"}
        missing = required - keys
        if missing:
            raise ValueError(f"DeepPBS NPZ missing keys: {sorted(missing)}")
        pwm = np.asarray(data["P"], dtype=float)
        seq_array = np.asarray(data["Seq"])

    if pwm.ndim != 2 or pwm.shape[1] != 4:
        raise ValueError(f"P must have shape (L, 4), got {pwm.shape}")
    if not np.isfinite(pwm).all() or (pwm < 0).any():
        raise ValueError("P contains non-finite or negative values")
    if not np.allclose(pwm.sum(axis=1), 1.0, atol=1e-5):
        raise ValueError("P rows are not normalized probabilities")
    sequence = _one_hot_to_sequence(seq_array)
    if len(sequence) != pwm.shape[0]:
        raise ValueError(f"P/Seq length mismatch: {pwm.shape[0]} vs {len(sequence)}")
    return {
        "pwm": pwm,
        "sequence": sequence,
        "length": int(pwm.shape[0]),
        "base_order": DEEPPBS_BASE_ORDER,
        "source_npz": str(path),
    }


def pwm_oriented_score(
    pwm: np.ndarray,
    sequence: str,
    *,
    epsilon: float = PRIMARY_EPSILON,
) -> tuple[float, int]:
    """Score a 7-mer against all contiguous DeepPBS PWM windows.

    Returns ``(best_log_probability, best_zero_based_offset)``.  This is a
    sequence-ranking proxy derived from DeepPBS probabilities, not affinity.
    """

    pwm = np.asarray(pwm, dtype=float)
    sequence = sequence.upper()
    if pwm.ndim != 2 or pwm.shape[1] != 4:
        raise ValueError(f"PWM must have shape (L, 4), got {pwm.shape}")
    if len(sequence) != 7 or any(base not in BASE_TO_INDEX for base in sequence):
        raise ValueError(f"Expected an ACGT 7-mer, got {sequence!r}")
    if pwm.shape[0] < 7:
        raise ValueError(f"DeepPBS PWM length must be at least 7, got {pwm.shape[0]}")
    indices = np.array([BASE_TO_INDEX[base] for base in sequence], dtype=int)
    log_pwm = np.log(np.clip(pwm, epsilon, 1.0))
    scores = [float(log_pwm[offset + np.arange(7), indices].sum()) for offset in range(pwm.shape[0] - 6)]
    best_offset = int(np.argmax(scores))
    return scores[best_offset], best_offset


def pwm_rc_class_score(
    pwm: np.ndarray,
    canonical_sequence: str,
    *,
    epsilon: float = PRIMARY_EPSILON,
) -> dict[str, object]:
    """Score one canonical RC class using the fixed max-over-orientations rule."""

    canonical_sequence = canonical_sequence.upper()
    if canonical_rc(canonical_sequence) != canonical_sequence:
        raise ValueError(f"Input is not canonical RC sequence: {canonical_sequence!r}")
    rc = reverse_complement(canonical_sequence)
    forward_score, forward_offset = pwm_oriented_score(pwm, canonical_sequence, epsilon=epsilon)
    reverse_score, reverse_offset = pwm_oriented_score(pwm, rc, epsilon=epsilon)
    if reverse_score > forward_score:
        best_score, best_orientation, best_offset = reverse_score, "reverse_complement", reverse_offset
    else:
        best_score, best_orientation, best_offset = forward_score, "forward", forward_offset
    return {
        "canonical_7mer": canonical_sequence,
        "oriented_7mer": rc if best_orientation == "reverse_complement" else canonical_sequence,
        "prediction_score": float(best_score),
        "prediction_orientation": best_orientation,
        "best_offset": int(best_offset),
        "forward_score": float(forward_score),
        "reverse_complement_score": float(reverse_score),
    }
