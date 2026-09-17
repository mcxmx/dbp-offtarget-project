from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .audit import rc_representation_audit, synthetic_conditionality_benchmark


@dataclass(frozen=True)
class GateResult:
    name: str
    passed: bool
    observed: str
    threshold: str


def evaluate_pre_replay_gates() -> pd.DataFrame:
    rc = rc_representation_audit().iloc[0]
    synthetic = synthetic_conditionality_benchmark()
    gates = [
        GateResult(
            "rc_class_count",
            int(rc["n_canonical_rc_classes"]) == 8192,
            str(int(rc["n_canonical_rc_classes"])),
            "8192",
        ),
        GateResult(
            "v06_rc_class_injectivity",
            bool(rc["v06_rc_class_injective"]),
            f"collisions={int(rc['v06_cross_rc_class_representation_collisions'])}",
            "zero cross-class collisions",
        ),
        GateResult(
            "synthetic_all_rankings_correct",
            bool(synthetic["ranking_correct"].all()),
            f"{int(synthetic['ranking_correct'].sum())}/{len(synthetic)}",
            "4/4",
        ),
        GateResult(
            "synthetic_protein_conditionality",
            bool(synthetic["protein_shuffle_not_invariant"].all()),
            f"min corr={synthetic['protein_shuffle_correlation'].min():.6f}",
            "positive correlation < 0.995",
        ),
        GateResult(
            "synthetic_target_conditionality",
            bool(synthetic["target_shuffle_not_invariant"].all()),
            f"min corr={synthetic['target_shuffle_correlation'].min():.6f}",
            "positive correlation < 0.995",
        ),
    ]
    return pd.DataFrame([gate.__dict__ for gate in gates])


def assert_pre_replay_gates() -> pd.DataFrame:
    gates = evaluate_pre_replay_gates()
    if not gates["passed"].all():
        failed = gates.loc[~gates["passed"], "name"].tolist()
        raise RuntimeError(f"v0.6 NO-GO before real-data replay: {failed}")
    return gates


def summarize_replay_gates(
    primary: pd.DataFrame,
    protein_shuffle: pd.DataFrame,
    target_shuffle: pd.DataFrame,
    health: pd.DataFrame,
) -> pd.DataFrame:
    model_rows = primary.loc[primary["model"].isin(("M0", "M1", "M1c", "M2", "M3"))]
    rows = [
        GateResult(
            "held_out_spearman_available",
            int(model_rows["spearman"].notna().sum()) == 105,
            f"{int(model_rows['spearman'].notna().sum())}/105",
            "105 complete protein-seed rows",
        ),
        GateResult(
            "protein_shuffle_no_collapse",
            float(protein_shuffle["prediction_correlation"].median()) <= 0.995,
            f"median corr={protein_shuffle['prediction_correlation'].median():.6f}",
            "median prediction correlation <= 0.995",
        ),
        GateResult(
            "target_shuffle_no_collapse",
            float(target_shuffle["prediction_correlation"].median()) <= 0.995,
            f"median corr={target_shuffle['prediction_correlation'].median():.6f}",
            "median prediction correlation <= 0.995",
        ),
        GateResult(
            "prediction_variance_positive",
            bool((health["prediction_variance"] > 0).all()),
            f"min variance={health['prediction_variance'].min():.6g}",
            "> 0 for every run",
        ),
        GateResult(
            "no_nan_inf",
            bool((health["nan_inf_count"] == 0).all()),
            f"max NaN/Inf={int(health['nan_inf_count'].max())}",
            "0",
        ),
    ]
    return pd.DataFrame([gate.__dict__ for gate in rows])
