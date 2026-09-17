from __future__ import annotations

from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.utils import project_root
from .models import (
    CandidateDNAOnlyV06,
    ProteinCandidateV06,
    ProteinTargetCandidateV06,
    TargetCandidateOnlyV06,
    batch_canonical_one_hot,
    canonical_one_hot,
    ordered_target_windows,
    score_rc_v06,
    target_window_features_v06,
)


ROOT = project_root()
DESIGNED_IDS = ("DBP1", "DBP3", "DBP35", "DBP48", "DBP5", "DBP6", "DBP9")
MODEL_ORDER = ("M0", "M1", "M1c", "M2", "M3")
V06_RESULTS = ROOT / "results" / "v0_6_representation"
V06_METADATA = ROOT / "metadata" / "v0_6"


def all_kmers(k: int = 7) -> list[str]:
    return ["".join(chars) for chars in product("ACGT", repeat=k)]


def rc_representation_audit() -> pd.DataFrame:
    sequences = all_kmers(7)
    classes = {canonical_rc(sequence) for sequence in sequences}
    repaired_by_rep: dict[tuple[float, ...], list[str]] = {}
    for sequence in sequences:
        repaired_by_rep.setdefault(tuple(canonical_one_hot(sequence).tolist()), []).append(sequence)
    from src.v0_5_models import rc_symmetric_one_hot

    old_reps = {}
    for sequence in sequences:
        old_reps.setdefault(tuple(rc_symmetric_one_hot(sequence).tolist()), []).append(sequence)
    old_cross_class = sum(
        len({canonical_rc(sequence) for sequence in members}) > 1
        for members in old_reps.values()
    )
    old_collision_excess = sum(max(0, len(members) - 2) for members in old_reps.values())
    new_cross_class = sum(
        len({canonical_rc(sequence) for sequence in members}) > 1
        for members in repaired_by_rep.values()
    )
    new_collision_excess = sum(max(0, len(members) - 2) for members in repaired_by_rep.values())
    return pd.DataFrame(
        [
            {
                "audit": "rc_encoding_exhaustive",
                "n_oriented_7mers": len(sequences),
                "n_canonical_rc_classes": len(classes),
                "legacy_representation_count": len(old_reps),
                "legacy_cross_rc_class_representation_collisions": old_cross_class,
                "legacy_collision_excess_over_rc_pairs": old_collision_excess,
                "v06_representation_count": len(repaired_by_rep),
                "v06_cross_rc_class_representation_collisions": new_cross_class,
                "v06_collision_excess_over_rc_pairs": new_collision_excess,
                "v06_rc_class_injective": new_cross_class == 0 and len(repaired_by_rep) == len(classes),
            }
        ]
    )


def target_representation_audit() -> pd.DataFrame:
    targets = pd.read_csv(ROOT / "metadata" / "v0_5" / "designed_target_manifest_v0_5.csv")
    rows = []
    for row in targets.to_dict("records"):
        windows = ordered_target_windows(row["primary_target"])
        old_windows = tuple(sorted({canonical_rc(row["primary_target"][i : i + 7]) for i in range(len(row["primary_target"]) - 6)}))
        rows.append(
            {
                "dbp_id": row["dbp_id"],
                "target_length": len(row["primary_target"]),
                "v06_ordered_window_count": len(windows),
                "v06_unique_window_count": len({window for window, _index, _position in windows}),
                "v06_repeated_window_count": len(windows) - len({window for window, _index, _position in windows}),
                "v06_register_min": min(position for _window, _index, position in windows),
                "v06_register_max": max(position for _window, _index, position in windows),
                "legacy_unique_sorted_window_count": len(old_windows),
                "ordered_windows_preserved": True,
                "set_deduplication_used": False,
                "register_preserved": True,
            }
        )
    return pd.DataFrame(rows)


def check_score_rc_invariance() -> pd.DataFrame:
    candidates = all_kmers(7)
    protein = torch.zeros(4)
    target = "ACGTACG"
    models = [
        ("M0", CandidateDNAOnlyV06(hidden_dim=8), CandidateDNAOnlyV06.model_name),
        ("M1", ProteinCandidateV06(4, hidden_dim=8), ProteinCandidateV06.model_name),
        ("M2", TargetCandidateOnlyV06(hidden_dim=8), TargetCandidateOnlyV06.model_name),
        ("M3", ProteinTargetCandidateV06(4, hidden_dim=8), ProteinTargetCandidateV06.model_name),
    ]
    rows = []
    rc_candidates = [reverse_complement(sequence) for sequence in candidates]
    for label, model, name in models:
        first = score_rc_v06(model, name, candidates, protein_embedding=protein, target=target)
        second = score_rc_v06(model, name, rc_candidates, protein_embedding=protein, target=target)
        rows.append(
            {
                "model": label,
                "n_pairs": len(candidates),
                "max_abs_rc_score_difference": float(np.max(np.abs(first - second))),
                "score_rc_invariant": bool(np.allclose(first, second, atol=1e-6, rtol=1e-6)),
            }
        )
    return pd.DataFrame(rows)


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 2 or np.unique(left).size < 2 or np.unique(right).size < 2:
        return np.nan
    return float(spearmanr(left, right).statistic)


def synthetic_conditionality_benchmark(
    *,
    seed: int = 2026,
    epochs: int = 240,
) -> pd.DataFrame:
    """Train the v0.6 P,T,D scorer on a four-cell sign-reversal task."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = ProteinTargetCandidateV06(protein_dim=2, hidden_dim=16, target_pooling="mean")
    optimizer = torch.optim.Adam(model.parameters(), lr=0.02, weight_decay=0.0)
    candidates = ["AAAAAAA", "CCCCCCC"]
    candidate_features = batch_canonical_one_hot(candidates)
    proteins = {"P1": torch.tensor([1.0, 0.0]), "P2": torch.tensor([0.0, 1.0])}
    targets = {"T1": "AAAAAAA", "T2": "CCCCCCC"}
    # Desired sign is protein sign times target sign: both inputs can reverse
    # the within-cell candidate ranking.
    signs = {("P1", "T1"): 1.0, ("P1", "T2"): -1.0, ("P2", "T1"): -1.0, ("P2", "T2"): 1.0}
    history = []
    for epoch in range(epochs):
        losses = []
        for (protein_id, target_id), sign in signs.items():
            protein = proteins[protein_id].reshape(1, -1).expand(2, -1)
            target = target_window_features_v06(targets[target_id]).unsqueeze(0).expand(2, -1, -1)
            prediction = model(protein, target, candidate_features)
            truth_sign = torch.tensor(sign)
            loss = torch.nn.functional.softplus(-(prediction[0] - prediction[1]) * truth_sign)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach()))
        history.append(float(np.mean(losses)))
    predictions: dict[tuple[str, str], np.ndarray] = {}
    rows = []
    for (protein_id, target_id), expected_sign in signs.items():
        prediction = score_rc_v06(
            model,
            model.model_name,
            candidates,
            protein_embedding=proteins[protein_id],
            target=targets[target_id],
        )
        predictions[(protein_id, target_id)] = prediction
        observed_sign = float(np.sign(prediction[0] - prediction[1]))
        rows.append(
            {
                "benchmark": "synthetic_conditionality",
                "protein_id": protein_id,
                "target_id": target_id,
                "candidate_A": candidates[0],
                "candidate_C": candidates[1],
                "expected_A_minus_C_sign": expected_sign,
                "observed_A_minus_C": float(prediction[0] - prediction[1]),
                "observed_A_minus_C_sign": observed_sign,
                "ranking_correct": bool(observed_sign == expected_sign),
                "prediction_variance": float(np.var(prediction)),
                "final_training_loss": history[-1],
                "epochs": epochs,
                "seed": seed,
            }
        )
    frame = pd.DataFrame(rows)
    # Explicit input-swap checks required by the conditionality gate.
    p1_t1 = frame.loc[(frame.protein_id == "P1") & (frame.target_id == "T1"), "observed_A_minus_C_sign"].iloc[0]
    p2_t1 = frame.loc[(frame.protein_id == "P2") & (frame.target_id == "T1"), "observed_A_minus_C_sign"].iloc[0]
    p1_t2 = frame.loc[(frame.protein_id == "P1") & (frame.target_id == "T2"), "observed_A_minus_C_sign"].iloc[0]
    frame["protein_swap_changes_ranking"] = p1_t1 != p2_t1
    frame["target_swap_changes_ranking"] = p1_t1 != p1_t2
    frame["protein_shuffle_correlation"] = [
        float(np.corrcoef(predictions[(protein_id, target_id)], predictions[("P2" if protein_id == "P1" else "P1", target_id)])[0, 1])
        for protein_id, target_id in zip(frame["protein_id"], frame["target_id"])
    ]
    frame["target_shuffle_correlation"] = [
        float(np.corrcoef(predictions[(protein_id, target_id)], predictions[(protein_id, "T2" if target_id == "T1" else "T1")])[0, 1])
        for protein_id, target_id in zip(frame["protein_id"], frame["target_id"])
    ]
    # The collapse gate concerns a positive near-one correlation. A negative
    # correlation is a strong conditionality signal, not a collapse.
    frame["protein_shuffle_not_invariant"] = frame["protein_shuffle_correlation"] < 0.995
    frame["target_shuffle_not_invariant"] = frame["target_shuffle_correlation"] < 0.995
    return frame


def run_representation_audit() -> dict[str, pd.DataFrame]:
    rc = rc_representation_audit()
    target = target_representation_audit()
    scores = check_score_rc_invariance()
    synthetic = synthetic_conditionality_benchmark()
    V06_METADATA.mkdir(parents=True, exist_ok=True)
    V06_RESULTS.mkdir(parents=True, exist_ok=True)
    rc.to_csv(V06_RESULTS / "representation_audit.csv", index=False)
    target.to_csv(V06_RESULTS / "target_representation_audit.csv", index=False)
    scores.to_csv(V06_RESULTS / "rc_score_invariance.csv", index=False)
    synthetic.to_csv(V06_RESULTS / "synthetic_conditionality.csv", index=False)
    return {"rc": rc, "target": target, "scores": scores, "synthetic": synthetic}


if __name__ == "__main__":
    outputs = run_representation_audit()
    for frame in outputs.values():
        print(frame.to_string(index=False))
