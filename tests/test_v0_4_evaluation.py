from pathlib import Path

import numpy as np
import pandas as pd

from src.v0_4_evaluation import compute_ranking_metrics


ROOT = Path(__file__).resolve().parents[1]


def test_ranking_metrics_use_high_score_as_better():
    df = pd.DataFrame(
        {
            "canonical_7mer": ["AAAAAAA", "AAAAAAC", "AAAAAAG", "AAAAAAT"],
            "truth": [0.1, 0.2, 0.3, 0.4],
            "perfect": [1.0, 2.0, 3.0, 4.0],
            "reversed": [4.0, 3.0, 2.0, 1.0],
        }
    )
    perfect = compute_ranking_metrics(df, "truth", "perfect")
    reversed_metric = compute_ranking_metrics(df, "truth", "reversed")
    assert np.isclose(perfect.spearman, 1.0)
    assert np.isclose(reversed_metric.spearman, -1.0)
    assert perfect.top1pct_recovery == 1.0
    assert reversed_metric.top1pct_recovery == 0.0


def test_v0_4_metrics_are_per_protein_then_macro_summarized():
    per = pd.read_csv(ROOT / "results" / "v0_4" / "tables" / "baseline_performance_per_protein.csv")
    macro = pd.read_csv(ROOT / "results" / "v0_4" / "tables" / "baseline_performance_macro.csv")
    assert set(per["protein_id"]) == {"DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"}
    assert per.groupby("baseline")["protein_id"].nunique().eq(7).all()
    assert macro["n_proteins"].max() <= 7
    sequence_kmer3 = macro[(macro["baseline"] == "sequence_kmer3") & (macro["metric"] == "spearman")].iloc[0]
    assert sequence_kmer3["n_proteins"] == 7
    assert 0.0 < sequence_kmer3["median"] < 0.5


def test_missing_external_predictions_are_not_filled_with_zero():
    scored = pd.read_parquet(ROOT / "data" / "processed" / "v0_4" / "v0_4_scored_candidates.parquet")
    assert scored["deeppbs_score"].isna().all()
    assert scored["simple_pc_score"].isna().all()
    non_evaluable = scored[~scored["protein_id"].isin(["DBP35", "DBP48"])]
    assert non_evaluable["structural_ppm_score"].isna().all()
