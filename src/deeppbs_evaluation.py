"""Project-side DeepPBS output parsing, scoring, and ranking evaluation."""

from __future__ import annotations

from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

from src.deeppbs_adapter import parse_prediction_npz, pwm_rc_class_score
from src.sequence_equivalence import canonical_rc
from src.v0_4_evaluation import (
    bootstrap_metric_ci,
    compute_ranking_metrics,
)


DNA_ALPHABET = "ACGT"
DESIGNED_PROTEINS = ("DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48")


def canonical_7mer_universe() -> list[str]:
    """Return the deterministic 8,192-element RC-class universe."""

    return sorted(
        {
            canonical_rc("".join(bases))
            for bases in product(DNA_ALPHABET, repeat=7)
        }
    )


def load_designed_experimental_units(project_root: str | Path) -> pd.DataFrame:
    """Load the frozen v0.3.1 designed-DBP RC-class PBM table."""

    path = (
        Path(project_root)
        / "data"
        / "processed"
        / "v0_3_1"
        / "designed_dbp_upbm_rc_class_v0_3_1.parquet"
    )
    df = pd.read_parquet(path).rename(
        columns={"experimental_escore_consensus": "experimental_score"}
    )
    required = {"protein_id", "canonical_7mer", "experimental_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Experimental table missing columns: {sorted(missing)}")
    if len(df) != len(DESIGNED_PROTEINS) * len(canonical_7mer_universe()):
        raise ValueError(f"Unexpected experimental unit count: {len(df)}")
    return df


def score_prediction_npz(
    npz_path: str | Path,
    protein_id: str,
    *,
    structure_id: str,
    structure_source: str,
    model_version: str,
    weight_id: str,
    overlap_status: str,
) -> pd.DataFrame:
    """Map one official DeepPBS prediction NPZ to all canonical 7-mers."""

    parsed = parse_prediction_npz(npz_path)
    rows = []
    for canonical_sequence in canonical_7mer_universe():
        score = pwm_rc_class_score(parsed["pwm"], canonical_sequence)
        rows.append(
            {
                "protein_id": protein_id,
                "canonical_7mer": canonical_sequence,
                "deeppbs_score": score["prediction_score"],
                "deeppbs_forward_score": score["forward_score"],
                "deeppbs_reverse_complement_score": score[
                    "reverse_complement_score"
                ],
                "deeppbs_orientation": score["prediction_orientation"],
                "deeppbs_best_offset": score["best_offset"],
                "structure_id": structure_id,
                "structure_source": structure_source,
                "structure_sequence": parsed["sequence"],
                "structure_length": parsed["length"],
                "model_version": model_version,
                "weight_id": weight_id,
                "overlap_status": overlap_status,
                "prediction_type": "DeepPBS PWM-derived log-probability proxy",
                "scoring_protocol": (
                    "max over contiguous 7-mer windows and candidate orientations; "
                    "sum log(P+1e-9)"
                ),
            }
        )
    return pd.DataFrame(rows)


def evaluate_per_protein(
    experimental: pd.DataFrame,
    predictions: pd.DataFrame,
) -> pd.DataFrame:
    """Compute ranking metrics separately for each protein."""

    rows = []
    for protein_id in DESIGNED_PROTEINS:
        exp = experimental[experimental["protein_id"] == protein_id][
            ["protein_id", "canonical_7mer", "experimental_score"]
        ]
        pred = predictions[predictions["protein_id"] == protein_id][
            ["protein_id", "canonical_7mer", "deeppbs_score"]
        ]
        merged = exp.merge(pred, on=["protein_id", "canonical_7mer"], how="left")
        metrics = compute_ranking_metrics(
            merged,
            truth_col="experimental_score",
            prediction_col="deeppbs_score",
        )
        rows.append(
            {
                "protein_id": protein_id,
                "baseline": "DeepPBS",
                "status": (
                    "evaluated"
                    if merged["deeppbs_score"].notna().any()
                    else "not_evaluable_missing_prediction"
                ),
                "spearman": metrics.spearman,
                "ndcg_1pct": metrics.ndcg_1pct,
                "ndcg_5pct": metrics.ndcg_5pct,
                "pairwise_accuracy": metrics.pairwise_accuracy,
                "top1pct_recovery": metrics.top1pct_recovery,
                "n_rc_classes": metrics.n_rc_classes,
                "n_expected_rc_classes": len(canonical_7mer_universe()),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_per_protein(
    experimental: pd.DataFrame,
    predictions: pd.DataFrame,
    *,
    n_bootstrap: int = 500,
    seed: int = 42,
) -> pd.DataFrame:
    """Bootstrap ranking metrics using canonical RC classes as units."""

    rows = []
    for protein_id in DESIGNED_PROTEINS:
        merged = experimental[experimental["protein_id"] == protein_id][
            ["protein_id", "canonical_7mer", "experimental_score"]
        ].merge(
            predictions[predictions["protein_id"] == protein_id][
                ["protein_id", "canonical_7mer", "deeppbs_score"]
            ],
            on=["protein_id", "canonical_7mer"],
            how="left",
        )
        if merged["deeppbs_score"].notna().sum() == 0:
            rows.append(
                {
                    "protein_id": protein_id,
                    "baseline": "DeepPBS",
                    "metric": "spearman",
                    "n_units": 0,
                    "bootstrap_mean": np.nan,
                    "ci_2.5": np.nan,
                    "ci_97.5": np.nan,
                    "bootstrap_unit": "canonical_rc_class",
                    "n_bootstrap": n_bootstrap,
                }
            )
            continue
        mean, lower, upper = bootstrap_metric_ci(
            merged,
            truth_col="experimental_score",
            prediction_col="deeppbs_score",
            metric="spearman",
            n_bootstrap=n_bootstrap,
            seed=seed + DESIGNED_PROTEINS.index(protein_id),
        )
        rows.append(
            {
                "protein_id": protein_id,
                "baseline": "DeepPBS",
                "metric": "spearman",
                "n_units": int(merged["deeppbs_score"].notna().sum()),
                "bootstrap_mean": mean,
                "ci_2.5": lower,
                "ci_97.5": upper,
                "bootstrap_unit": "canonical_rc_class",
                "n_bootstrap": n_bootstrap,
            }
        )
    return pd.DataFrame(rows)
