from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/v1_3"

def corr(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0: return np.nan
    return float(spearmanr(x, y).statistic)

def pairwise(x, y):
    vals = []
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            dx, dy = x[i] - x[j], y[i] - y[j]
            if dx == 0: continue
            vals.append(1.0 if dx * dy > 0 else 0.5 if dy == 0 else 0.0)
    return float(np.mean(vals)) if vals else np.nan

def variance_rows(df, source_type, cohort):
    rows = []
    for protein, g in df.groupby("protein_id", sort=True):
        y = g.value.to_numpy(float); mu = y.mean()
        alpha = g.groupby("position").value.transform("mean").to_numpy(float) - mu
        eps = y - mu - alpha; total = float(np.var(y)); between = float(np.var(alpha)); within = float(np.var(eps))
        rows.append({"protein": protein, "cohort": cohort, "source_type": source_type, "total_variance": total, "position_variance": between, "identity_variance": within, "position_fraction": between / total if total else np.nan, "identity_fraction": within / total if total else np.nan, "n_mutations": len(g), "n_positions": g.position.nunique()})
    return rows

def main():
    exp = pd.read_parquet(ROOT / "data/processed/v1_3_competition_mutations_official_mean.parquet")
    exp = exp[["protein_id", "position", "wt_base", "mutant_base", "normalized_effect"]].rename(columns={"normalized_effect": "value"})
    pred = pd.read_csv(ROOT / "results/v1_1/05_assay_aligned_eval/deeppbs_per_mutation_predictions.tsv", sep="	")
    pred = pred.rename(columns={"protein": "protein_id", "assay_position": "position", "mut_base": "mutant_base", "predicted_delta_logP": "value"})
    pred = pred[["protein_id", "position", "wt_base", "mutant_base", "value"]]
    pred = pred[pred.protein_id.isin(["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP035", "DBP048"])]
    hold = pd.read_csv(OUT / "holdout_predictions_revealed.tsv", sep="	")
    rows = variance_rows(exp[exp.protein_id.isin(pred.protein_id)], "experimental_competition", "development")
    rows += variance_rows(pred, "DeepPBS_prediction", "development")
    rows += variance_rows(exp[exp.protein_id.isin(["DBP023", "DBP056", "DBP062"])], "experimental_competition", "prospective_holdout")
    hp = hold[["protein_id", "position", "wt_base", "mutant_base", "predicted_effect"]].rename(columns={"predicted_effect": "value"})
    rows += variance_rows(hp, "DeepPBS_prediction", "prospective_holdout")
    pd.DataFrame(rows).to_csv(OUT / "variance_decomposition_per_protein.tsv", sep="	", index=False, na_rep="NA")
    cf = []
    joined = pred.merge(exp, on=["protein_id", "position", "wt_base", "mutant_base"], suffixes=("_pred", "_exp"))
    for cohort, d in [("development", joined), ("prospective_holdout", hold.rename(columns={"predicted_effect": "value_pred", "normalized_effect": "value_exp"}))]:
        for protein, g in d.groupby("protein_id", sort=True):
            g = g.sort_values(["position", "mutant_base"]); full = g.value_pred.to_numpy(float); ex = g.value_exp.to_numpy(float)
            pos = g.groupby("position").value_pred.transform("mean").to_numpy(float)
            er = g.groupby("position").value_exp.transform(lambda z: z - z.mean()).to_numpy(float)
            pr = g.groupby("position").value_pred.transform(lambda z: z - z.mean()).to_numpy(float)
            cf.append({"protein": protein, "cohort": cohort, "full_global_spearman": corr(ex, full), "position_only_global_spearman": corr(ex, pos), "delta_global": corr(ex, pos) - corr(ex, full), "full_identity_spearman": corr(er, pr), "full_pairwise_accuracy": float(np.nanmean([pairwise(a.value_exp.to_numpy(float), a.value_pred.to_numpy(float)) for _, a in g.groupby("position")]))})
    pd.DataFrame(cf).to_csv(OUT / "position_only_counterfactual.tsv", sep="	", index=False, na_rep="NA")
    pd.DataFrame([{"step":"positive_control", "status":"PASS_CACHED_SUCCESS", "evidence":"DBP035.npz and DBP035.npz_predict.npz exist", "route":"existing frozen local pipeline"}, {"step":"2STT_MODEL1_local", "status":"FAILED_ENVIRONMENT_MISSING_TORCH_CLUSTER", "evidence":"process_co_crystal.py exit 1; ModuleNotFoundError: torch_cluster; no NPZ", "route":"Route 1"}, {"step":"2STT_MODEL1_docker", "status":"NOT_ATTEMPTED_DOCKER_UNAVAILABLE_CACHED", "evidence":"existing audit records Docker exit 127 and WSL exit 50", "route":"Route 2"}, {"step":"2STT_final", "status":"NOT_EVALUABLE_PRIMARY_UNCHANGED", "evidence":"no rescue NPZ; labels revealed before rescue", "route":"bounded rescue stopped"}]).to_csv(OUT / "natural_tf_2stt_rescue_status.tsv", sep="	", index=False)

if __name__ == "__main__": main()
