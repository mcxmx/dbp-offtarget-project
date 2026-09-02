from __future__ import annotations

import itertools
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.sequence_equivalence import canonical_rc, reverse_complement


DNA_BASES = "ACGT"
DNA_TO_IDX = {base: idx for idx, base in enumerate(DNA_BASES)}
IDX_TO_DNA = {idx: base for base, idx in DNA_TO_IDX.items()}
DNA_RESIDUE_TO_BASE = {"DA": "A", "DC": "C", "DG": "G", "DT": "T"}


@dataclass(frozen=True)
class StructuralRun:
    protein_id: str
    structure_id: str
    chain_label: str
    start_resnum: int
    end_resnum: int
    sequence: str
    log_probs: np.ndarray  # shape [L, 4]


def _inverse_restype_map(restype_to_int: dict) -> dict[int, str]:
    inverse: dict[int, str] = {}
    for residue, idx in restype_to_int.items():
        if residue in DNA_RESIDUE_TO_BASE and idx not in inverse:
            inverse[idx] = DNA_RESIDUE_TO_BASE[residue]
    return inverse


def load_nampnn_npz(npz_path: str | Path) -> dict:
    data = np.load(npz_path, allow_pickle=True)
    restype_to_int = data["restype_to_int"].item()
    inverse = _inverse_restype_map(restype_to_int)
    residues = np.asarray(data["encoded_residues"], dtype=str)
    chain_labels = np.asarray(data["chain_labels"])
    dna_mask = np.asarray(data["dna_mask"]).astype(bool)
    true_sequence = np.asarray(data["true_sequence"])
    predicted_ppm = np.asarray(data["predicted_ppm"], dtype=float)
    dna_rows = []
    for idx, is_dna in enumerate(dna_mask):
        if not is_dna:
            continue
        match = re.match(r"^([A-Za-z]+)(\d+)$", residues[idx])
        if not match:
            continue
        chain_label, resnum_text = match.groups()
        resnum = int(resnum_text)
        base = inverse.get(int(true_sequence[idx]), "N")
        probs = predicted_ppm[idx, [restype_to_int["DA"], restype_to_int["DC"], restype_to_int["DG"], restype_to_int["DT"]]]
        dna_rows.append(
            {
                "index": idx,
                "chain_label": chain_label,
                "resnum": resnum,
                "base": base,
                "log_probs": np.log(np.clip(probs.astype(float), 1e-12, 1.0)),
                "chain_numeric_label": int(chain_labels[idx]),
            }
        )
    return {
        "residues": residues,
        "dna_rows": dna_rows,
        "restype_to_int": restype_to_int,
    }


def contiguous_runs(dna_rows: list[dict]) -> list[list[dict]]:
    runs: list[list[dict]] = []
    by_chain: dict[str, list[dict]] = {}
    for row in dna_rows:
        by_chain.setdefault(row["chain_label"], []).append(row)
    for chain_label, rows in by_chain.items():
        rows = sorted(rows, key=lambda r: r["resnum"])
        if not rows:
            continue
        current = [rows[0]]
        for row in rows[1:]:
            prev = current[-1]
            if row["resnum"] == prev["resnum"] + 1:
                current.append(row)
            else:
                runs.append(current)
                current = [row]
        runs.append(current)
    return runs


def score_sequence_against_runs(sequence: str, runs: list[list[dict]]) -> dict[str, object]:
    seq = sequence.upper()
    if len(seq) != 7 or any(base not in DNA_BASES for base in seq):
        raise ValueError(f"Expected ACGT 7-mer, got {sequence!r}")
    seq_idx = np.array([DNA_TO_IDX[base] for base in seq], dtype=int)
    rc_seq = reverse_complement(seq)
    rc_idx = np.array([DNA_TO_IDX[base] for base in rc_seq], dtype=int)

    def score_oriented(oriented_idx: np.ndarray) -> dict[str, object]:
        best = {
            "prediction_score": -np.inf,
            "best_chain_label": None,
            "best_window_start": None,
            "best_window_end": None,
            "best_window_sequence": None,
        }
        for run in runs:
            if len(run) < 7:
                continue
            log_probs = np.stack([row["log_probs"] for row in run], axis=0)
            chain_label = run[0]["chain_label"]
            for start in range(0, len(run) - 6):
                window = log_probs[start : start + 7, :]
                score = float(window[np.arange(7), oriented_idx].sum())
                if score > best["prediction_score"]:
                    best.update(
                        {
                            "prediction_score": score,
                            "best_chain_label": chain_label,
                            "best_window_start": int(run[start]["resnum"]),
                            "best_window_end": int(run[start + 6]["resnum"]),
                            "best_window_sequence": "".join(row["base"] for row in run[start : start + 7]),
                        }
                    )
        return best

    forward = score_oriented(seq_idx)
    reverse = score_oriented(rc_idx)
    if reverse["prediction_score"] > forward["prediction_score"]:
        return {
            "canonical_7mer": canonical_rc(seq),
            "oriented_7mer": rc_seq,
            "reverse_complement_7mer": seq,
            "prediction_score": reverse["prediction_score"],
            "best_chain_label": reverse["best_chain_label"],
            "best_window_start": reverse["best_window_start"],
            "best_window_end": reverse["best_window_end"],
            "best_window_sequence": reverse["best_window_sequence"],
            "prediction_orientation": "reverse_complement",
        }
    return {
        "canonical_7mer": canonical_rc(seq),
        "oriented_7mer": seq,
        "reverse_complement_7mer": rc_seq,
        "prediction_score": forward["prediction_score"],
        "best_chain_label": forward["best_chain_label"],
        "best_window_start": forward["best_window_start"],
        "best_window_end": forward["best_window_end"],
        "best_window_sequence": forward["best_window_sequence"],
        "prediction_orientation": "forward",
    }


def all_canonical_7mers() -> list[str]:
    canonical = sorted({canonical_rc("".join(chars)) for chars in itertools.product(DNA_BASES, repeat=7)})
    return canonical


def score_structural_ppm(
    npz_path: str | Path,
    protein_id: str,
    structure_id: str,
    model_version: str,
    prediction_type: str,
) -> pd.DataFrame:
    loaded = load_nampnn_npz(npz_path)
    runs = contiguous_runs(loaded["dna_rows"])
    canonical_7mers = all_canonical_7mers()
    rows = []
    for canonical in canonical_7mers:
        scored = score_sequence_against_runs(canonical, runs)
        rows.append(
            {
                "protein_id": protein_id,
                "canonical_7mer": scored["canonical_7mer"],
                "oriented_7mer": scored["oriented_7mer"],
                "reverse_complement_7mer": scored["reverse_complement_7mer"],
                "prediction_score": float(scored["prediction_score"]),
                "prediction_type": prediction_type,
                "structure_id": structure_id,
                "model_version": model_version,
                "best_chain_label": scored["best_chain_label"],
                "best_window_start": scored["best_window_start"],
                "best_window_end": scored["best_window_end"],
                "best_window_sequence": scored["best_window_sequence"],
                "prediction_orientation": scored["prediction_orientation"],
                "source_npz": str(Path(npz_path)),
            }
        )
    return pd.DataFrame(rows)

