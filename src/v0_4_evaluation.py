from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


@dataclass(frozen=True)
class RankingMetrics:
    spearman: float
    ndcg_1pct: float
    ndcg_5pct: float
    pairwise_accuracy: float
    top1pct_recovery: float
    n_rc_classes: int


def _rank_desc(values: pd.Series) -> pd.Series:
    return values.rank(method="average", ascending=False)


def ndcg_at_fraction(y_true: np.ndarray, y_score: np.ndarray, fraction: float) -> float:
    if len(y_true) == 0:
        return np.nan
    k = max(1, int(np.ceil(len(y_true) * fraction)))
    order = np.argsort(-y_score, kind="mergesort")[:k]
    ideal = np.argsort(-y_true, kind="mergesort")[:k]
    # E-scores can be negative, so use rank-derived non-negative gains rather
    # than treating the raw score as a relevance probability.
    true_order = pd.Series(y_true).rank(method="average", ascending=True).to_numpy()
    gains = true_order / len(true_order)
    discounts = 1.0 / np.log2(np.arange(2, k + 2))
    dcg = float(np.sum(gains[order] * discounts))
    idcg = float(np.sum(gains[ideal] * discounts))
    return dcg / idcg if idcg > 0 else np.nan


def top_fraction_recovery(y_true: np.ndarray, y_score: np.ndarray, fraction: float) -> float:
    if len(y_true) == 0:
        return np.nan
    k = max(1, int(np.ceil(len(y_true) * fraction)))
    true_top = set(np.argsort(-y_true, kind="mergesort")[:k])
    pred_top = set(np.argsort(-y_score, kind="mergesort")[:k])
    return len(true_top & pred_top) / k


def pairwise_ranking_accuracy(y_true: np.ndarray, y_score: np.ndarray, max_pairs: int = 20_000, seed: int = 42) -> float:
    n = len(y_true)
    if n < 2:
        return np.nan
    rng = np.random.default_rng(seed)
    total_pairs = n * (n - 1) // 2
    if total_pairs <= max_pairs:
        i, j = np.triu_indices(n, k=1)
    else:
        i = rng.integers(0, n, size=max_pairs)
        j = rng.integers(0, n - 1, size=max_pairs)
        j = np.where(j >= i, j + 1, j)
    true_delta = y_true[i] - y_true[j]
    score_delta = y_score[i] - y_score[j]
    valid = true_delta != 0
    if not np.any(valid):
        return np.nan
    product = true_delta[valid] * score_delta[valid]
    correct = (product > 0).sum() + 0.5 * (score_delta[valid] == 0).sum()
    return float(correct / valid.sum())


def compute_ranking_metrics(
    df: pd.DataFrame,
    truth_col: str,
    prediction_col: str,
    seed: int = 42,
) -> RankingMetrics:
    sub = df[[truth_col, prediction_col]].dropna()
    if sub.empty:
        return RankingMetrics(np.nan, np.nan, np.nan, np.nan, np.nan, 0)
    y_true = sub[truth_col].to_numpy(dtype=float)
    y_score = sub[prediction_col].to_numpy(dtype=float)
    if np.unique(y_score).size < 2 or np.unique(y_true).size < 2:
        spearman = np.nan
    else:
        spearman = float(spearmanr(y_true, y_score).statistic)
    return RankingMetrics(
        spearman=spearman,
        ndcg_1pct=ndcg_at_fraction(y_true, y_score, 0.01),
        ndcg_5pct=ndcg_at_fraction(y_true, y_score, 0.05),
        pairwise_accuracy=pairwise_ranking_accuracy(y_true, y_score, seed=seed),
        top1pct_recovery=top_fraction_recovery(y_true, y_score, 0.01),
        n_rc_classes=int(len(sub)),
    )


def add_percentiles(df: pd.DataFrame, score_col: str, out_col: str) -> pd.DataFrame:
    df = df.copy()
    df[out_col] = df.groupby("protein_id", group_keys=False)[score_col].rank(pct=True, ascending=True)
    return df


def bootstrap_metric_ci(
    df: pd.DataFrame,
    truth_col: str,
    prediction_col: str,
    metric: str = "spearman",
    n_bootstrap: int = 200,
    seed: int = 42,
) -> tuple[float, float, float]:
    clean = df[[truth_col, prediction_col, "canonical_7mer"]].dropna()
    if clean.empty:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    values = []
    indices = np.arange(len(clean))
    for _ in range(n_bootstrap):
        sample_idx = rng.choice(indices, size=len(indices), replace=True)
        sample = clean.iloc[sample_idx]
        if metric == "spearman":
            if sample[truth_col].nunique() < 2 or sample[prediction_col].nunique() < 2:
                values.append(np.nan)
            else:
                values.append(float(spearmanr(sample[truth_col], sample[prediction_col]).statistic))
        else:
            result = compute_ranking_metrics(sample, truth_col, prediction_col, seed=int(rng.integers(0, 2**31 - 1)))
            values.append(getattr(result, metric))
    arr = np.array(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return np.nan, np.nan, np.nan
    return float(np.nanmean(arr)), float(np.nanpercentile(arr, 2.5)), float(np.nanpercentile(arr, 97.5))


def macro_summary(per_protein: pd.DataFrame, metric_cols: list[str]) -> pd.DataFrame:
    rows = []
    for baseline, group in per_protein.groupby("baseline", sort=True):
        for metric in metric_cols:
            vals = group[metric].dropna().to_numpy(dtype=float)
            rows.append(
                {
                    "baseline": baseline,
                    "metric": metric,
                    "n_proteins_with_metric": int(vals.size),
                    "mean": float(np.mean(vals)) if vals.size else np.nan,
                    "median": float(np.median(vals)) if vals.size else np.nan,
                    "std": float(np.std(vals, ddof=1)) if vals.size > 1 else np.nan,
                }
            )
    return pd.DataFrame(rows)
