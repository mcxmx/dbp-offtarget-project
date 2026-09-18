from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/v1_3"
REPORTS = ROOT / "reports/v1.3"
FIGURES = ROOT / "figures/v1_3"
PREDICTIONS = RESULTS / "holdout_predictions_unscored.tsv"
LABELS = ROOT / "data/processed/v1_3_competition_mutations.parquet"
MANIFEST = ROOT / "artifacts/cache/file_manifest.json"
HOLDOUTS = ["DBP023", "DBP056", "DBP062"]
SEED = 1301
N_PERM = 2000
N_BOOT = 2000


def safe_corr(x, y, kind="spearman") -> float:
    a, b = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return float((spearmanr if kind == "spearman" else pearsonr)(a, b).statistic)


def pairwise_accuracy(experimental, predicted) -> float:
    exp, pred = np.asarray(experimental, float), np.asarray(predicted, float)
    scores = []
    for i in range(len(exp)):
        for j in range(i + 1, len(exp)):
            de, dp = exp[i] - exp[j], pred[i] - pred[j]
            if de == 0:
                continue
            scores.append(1.0 if de * dp > 0 else 0.5 if dp == 0 else 0.0)
    return float(np.mean(scores)) if scores else np.nan


def check_freeze() -> dict:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    key = PREDICTIONS.relative_to(ROOT).as_posix()
    frozen = manifest["files"][key]
    stat = PREDICTIONS.stat()
    if stat.st_size != frozen["size"] or stat.st_mtime_ns != frozen["mtime_ns"]:
        raise RuntimeError("Frozen prediction path/size/mtime changed; refusing reveal without a new documented freeze")
    return frozen


def position_metrics(table: pd.DataFrame) -> dict:
    groups = [g for _, g in table.groupby("position", sort=True)]
    exp_pos = np.array([np.mean(np.abs(g.normalized_effect)) for g in groups])
    pred_pos = np.array([np.mean(np.abs(g.predicted_effect)) for g in groups])
    exp_res = np.concatenate([(g.normalized_effect - g.normalized_effect.mean()).to_numpy(float) for g in groups])
    pred_res = np.concatenate([(g.predicted_effect - g.predicted_effect.mean()).to_numpy(float) for g in groups])
    pairwise = float(np.nanmean([pairwise_accuracy(g.normalized_effect, g.predicted_effect) for g in groups]))
    return {
        "n_positions": len(groups),
        "position_spearman": safe_corr(exp_pos, pred_pos),
        "position_pearson": safe_corr(exp_pos, pred_pos, "pearson"),
        "within_residual_spearman": safe_corr(exp_res, pred_res),
        "within_residual_pearson": safe_corr(exp_res, pred_res, "pearson"),
        "within_pairwise_accuracy": pairwise,
        "within_kendall_style_concordance": 2 * pairwise - 1,
    }


def stratified_rank_ci(table: pd.DataFrame, rng: np.random.Generator) -> tuple[float, float]:
    groups = [g for _, g in table.groupby("position", sort=True)]
    values = []
    for _ in range(500):
        sampled = pd.concat([g.iloc[rng.integers(0, len(g), len(g))] for g in groups], ignore_index=True)
        values.append(safe_corr(sampled.normalized_effect, sampled.predicted_effect))
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def evaluate_holdout() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    predictions = pd.read_csv(PREDICTIONS, sep="\t")
    labels = pd.read_parquet(LABELS)
    labels = labels[labels.protein_id.isin(HOLDOUTS)].copy()
    joined = labels.merge(
        predictions,
        left_on=["protein_id", "position", "wt_base", "mutant_base"],
        right_on=["protein", "position", "wt_base", "mutant_base"],
        how="left",
        validate="one_to_one",
    )
    if joined.predicted_effect.isna().any() or len(joined) != 126:
        raise ValueError("Holdout key join is incomplete")
    joined["experimental_identity_residual"] = joined.normalized_effect - joined.groupby(["protein_id", "position"]).normalized_effect.transform("mean")
    joined["predicted_identity_residual"] = joined.predicted_effect - joined.groupby(["protein_id", "position"]).predicted_effect.transform("mean")
    joined.to_csv(RESULTS / "holdout_predictions_revealed.tsv", sep="\t", index=False)

    per_rows, null_summary, null_distribution = [], [], []
    for protein_index, protein in enumerate(HOLDOUTS):
        q = joined[joined.protein_id == protein].sort_values(["position", "mutant_base"]).reset_index(drop=True)
        ci_low, ci_high = stratified_rank_ci(q, np.random.default_rng(SEED + protein_index))
        metrics = position_metrics(q)
        per_rows.append({
            "row_type": "protein",
            "protein_id": protein,
            "method": "DeepPBS",
            "holdout_status": "LOCKED_HOLDOUT_REVEALED_AFTER_FREEZE",
            "n_mutations": len(q),
            "unmatched_structure_predictions": int(len(predictions[predictions.protein == protein]) - len(q)),
            "global_spearman": safe_corr(q.normalized_effect, q.predicted_effect),
            "global_pearson": safe_corr(q.normalized_effect, q.predicted_effect, "pearson"),
            "global_rank_ci_low": ci_low,
            "global_rank_ci_high": ci_high,
            **metrics,
        })

        rng = np.random.default_rng(SEED + protein_index)
        exp = q.normalized_effect.to_numpy(float)
        pred = q.predicted_effect.to_numpy(float)
        groups = [g.index.to_numpy() for _, g in q.groupby("position", sort=True)]
        exp_pos = np.array([np.mean(np.abs(exp[idx])) for idx in groups])
        pred_pos = np.array([np.mean(np.abs(pred[idx])) for idx in groups])
        observed = position_metrics(q)
        values = {"position_spearman": [], "within_residual_spearman": [], "within_pairwise_accuracy": []}
        for permutation in range(N_PERM):
            values["position_spearman"].append(safe_corr(exp_pos, rng.permutation(pred_pos)))
            shuffled = pred.copy()
            for idx in groups:
                shuffled[idx] = rng.permutation(shuffled[idx])
            exp_res = np.concatenate([exp[idx] - np.mean(exp[idx]) for idx in groups])
            pred_res = np.concatenate([shuffled[idx] - np.mean(shuffled[idx]) for idx in groups])
            values["within_residual_spearman"].append(safe_corr(exp_res, pred_res))
            values["within_pairwise_accuracy"].append(float(np.nanmean([pairwise_accuracy(exp[idx], shuffled[idx]) for idx in groups])))
        definitions = [
            ("position_label", "position_spearman", 0.0, "two_sided"),
            ("within_position_base_label", "within_residual_spearman", 0.0, "two_sided"),
            ("within_position_base_label", "within_pairwise_accuracy", 0.5, "upper"),
        ]
        for null_type, metric, chance, tail in definitions:
            array = np.asarray(values[metric], float)
            array = array[np.isfinite(array)]
            obs = observed[metric]
            count = np.sum(np.abs(array - chance) >= abs(obs - chance)) if tail == "two_sided" else np.sum(array >= obs)
            null_summary.append({
                "protein_id": protein, "method": "DeepPBS", "null_type": null_type, "metric": metric,
                "observed": obs, "chance_baseline": chance, "null_median": np.median(array),
                "null_q025": np.quantile(array, 0.025), "null_q975": np.quantile(array, 0.975),
                "empirical_p": (1 + count) / (1 + len(array)), "tail": tail, "n_permutations": len(array),
            })
            null_distribution.extend({
                "protein_id": protein, "method": "DeepPBS", "null_type": null_type,
                "metric": metric, "permutation": i, "null_value": value,
            } for i, value in enumerate(array))

    per = pd.DataFrame(per_rows)
    metric_cols = ["global_spearman", "global_pearson", "position_spearman", "position_pearson", "within_residual_spearman", "within_residual_pearson", "within_pairwise_accuracy", "within_kendall_style_concordance"]
    median = {col: float(per[col].median()) for col in metric_cols}
    median.update({
        "row_type": "protein_median", "protein_id": "HOLDOUT_MEDIAN", "method": "DeepPBS",
        "holdout_status": "LOCKED_HOLDOUT_THREE_PROTEIN_SUMMARY", "n_mutations": int(per.n_mutations.sum()),
        "n_positions": int(per.n_positions.sum()), "unmatched_structure_predictions": int(per.unmatched_structure_predictions.sum()),
        "global_rank_ci_low": np.nan, "global_rank_ci_high": np.nan,
    })
    output = pd.concat([per, pd.DataFrame([median])], ignore_index=True)
    output.to_csv(RESULTS / "holdout_decomposition.tsv", sep="\t", index=False, na_rep="NA")

    summary_rows = []
    rng = np.random.default_rng(SEED)
    for metric in ["global_spearman", "position_spearman", "within_residual_spearman", "within_pairwise_accuracy"]:
        x = per[metric].to_numpy(float)
        draws = np.array([np.median(rng.choice(x, len(x), replace=True)) for _ in range(N_BOOT)])
        summary_rows.append({
            "analysis_set": "locked_holdout_three", "method": "DeepPBS", "metric": metric,
            "n_proteins": len(x), "median": np.median(x), "mean": np.mean(x),
            "ci_low": np.quantile(draws, 0.025), "ci_high": np.quantile(draws, 0.975),
            "chance_baseline": 0.5 if metric == "within_pairwise_accuracy" else 0.0,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(RESULTS / "holdout_decomposition_summary.tsv", sep="\t", index=False)
    pd.DataFrame(null_summary).to_csv(RESULTS / "holdout_permutation_nulls.tsv", sep="\t", index=False)
    pd.DataFrame(null_distribution).to_csv(RESULTS / "holdout_permutation_distributions.tsv", sep="\t", index=False)
    return per, summary, joined


def figure2_revised(per: pd.DataFrame) -> None:
    existing = pd.read_csv(RESULTS / "figure2_data.tsv", sep="\t")
    holdout = pd.DataFrame({
        "method": "DeepPBS", "protein_id": per.protein_id,
        "global_spearman": per.global_spearman, "position_spearman": per.position_spearman,
        "within_position_residual_spearman": per.within_residual_spearman,
        "within_position_pairwise_accuracy": per.within_pairwise_accuracy,
        "within_kendall_style_concordance": per.within_kendall_style_concordance,
        "n_mutations": per.n_mutations, "n_positions": per.n_positions,
        "holdout_status": "LOCKED_HOLDOUT", "analysis_status": "LOCKED_HOLDOUT",
        "mapping_status": "PREDICTION_PROTOCOL_FROZEN", "coverage_scope": "assay_matched_after_structure_wide_freeze",
        "global_chance": 0.0, "position_chance": 0.0, "within_residual_chance": 0.0,
        "within_pairwise_chance": 0.5, "within_kendall_chance": 0.0,
    })
    source = pd.concat([existing, holdout], ignore_index=True)
    source.to_csv(RESULTS / "figure2_revised_data.tsv", sep="\t", index=False)

    deep = source[(source.method == "DeepPBS") & (source.analysis_status != "SENSITIVITY_ONLY")].copy()
    deep["cohort"] = np.where(deep.holdout_status == "LOCKED_HOLDOUT", "Locked holdout", "Development-exposed")
    cols = ["global_spearman", "position_spearman", "within_position_residual_spearman"]
    labels = ["Global", "Position", "Identity residual"]
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    colors = {"Development-exposed": "#667085", "Locked holdout": "#C43D3D"}
    for cohort, group in deep.groupby("cohort"):
        for _, row in group.iterrows():
            ax.plot(range(3), [row[c] for c in cols], color=colors[cohort], alpha=0.3, lw=1)
            ax.scatter(range(3), [row[c] for c in cols], color=colors[cohort], alpha=0.75, s=26)
        med = [group[c].median() for c in cols]
        ax.plot(range(3), med, color=colors[cohort], lw=3, marker="o", ms=7, label=f"{cohort} median (n={len(group)})")
    ax.axhline(0, color="black", lw=0.8, alpha=0.5)
    ax.set_xticks(range(3), labels)
    ax.set_ylabel("Per-protein Spearman rho")
    ax.set_title("DeepPBS local-effect decomposition")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#E4E7EC", lw=0.7)
    fig.tight_layout()
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"Figure2_revised_development_holdout.{suffix}", dpi=300)
    plt.close(fig)


def figure3_reproducibility(per: pd.DataFrame) -> None:
    development = pd.read_csv(RESULTS / "figure2_data.tsv", sep="\t")
    development = development[(development.method == "DeepPBS") & (development.analysis_status == "PRIMARY")]
    combined = pd.concat([
        development[["protein_id", "global_spearman", "position_spearman", "within_position_residual_spearman"]].assign(cohort="development_exposed"),
        per[["protein_id", "global_spearman", "position_spearman", "within_residual_spearman"]].rename(columns={"within_residual_spearman": "within_position_residual_spearman"}).assign(cohort="locked_holdout"),
    ], ignore_index=True)
    rows = []
    metrics = {"GLOBAL": "global_spearman", "POSITION": "position_spearman", "IDENTITY": "within_position_residual_spearman"}
    for _, row in combined.iterrows():
        for level, metric in metrics.items():
            rows.append({"protein_id": row.protein_id, "cohort": row.cohort, "comparison": "DeepPBS_vs_competition", "level": level, "correlation": row[metric], "status": "EVALUABLE"})
    for protein in ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP023", "DBP035", "DBP048", "DBP056", "DBP062"]:
        for level in metrics:
            rows.append({"protein_id": protein, "cohort": "experimental_replicate", "comparison": "competition_R1_vs_R2", "level": level, "correlation": np.nan, "status": "NOT_EVALUABLE_PUBLIC_SOURCE_MEAN_ONLY"})
    source = pd.DataFrame(rows)
    source.to_csv(RESULTS / "figure3_reproducibility_data.tsv", sep="\t", index=False, na_rep="NA")

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 4.2), sharey=True)
    colors = {"development_exposed": "#667085", "locked_holdout": "#C43D3D"}
    rng = np.random.default_rng(SEED)
    for ax, level in zip(axes, metrics):
        evaluable = source[(source.level == level) & (source.comparison == "DeepPBS_vs_competition")]
        for cohort, group in evaluable.groupby("cohort"):
            x = 0 + rng.uniform(-0.08, 0.08, len(group))
            ax.scatter(x, group.correlation, color=colors[cohort], s=32, alpha=0.8, label=cohort.replace("_", " "))
        ax.axvspan(0.65, 1.35, color="#F2F4F7")
        ax.text(1, 0, "not evaluable\nmean-only source", ha="center", va="center", color="#667085", fontsize=9)
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)
        ax.set_xlim(-0.35, 1.35)
        ax.set_xticks([0, 1], ["DeepPBS vs\ncompetition", "Competition\nR1 vs R2"])
        ax.set_title(level.title())
        ax.grid(axis="y", color="#E4E7EC", lw=0.7)
    axes[0].set_ylabel("Spearman rho")
    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.91), ncol=2, frameon=False)
    fig.suptitle("Experimental reproducibility ceiling is absent from public source data", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.80))
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"Figure3_reproducibility_ceiling.{suffix}", dpi=300)
    plt.close(fig)


def figure4_dbp35opt() -> None:
    mutation = pd.read_csv(RESULTS / "dbp35opt_mutation_shift.tsv", sep="\t")
    position = pd.read_csv(RESULTS / "dbp35opt_position_shift.tsv", sep="\t")
    source = mutation.copy()
    source.to_csv(RESULTS / "figure4_dbp35opt_data.tsv", sep="\t", index=False)
    order = list("ACGT")

    def matrix(column):
        return mutation.pivot(index="mutant_base", columns="position", values=column).reindex(order).to_numpy(float)

    base = matrix("normalized_effect_dbp35")
    opt = matrix("normalized_effect_dbp35opt")
    delta = matrix("experimental_effect_shift")
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.4), gridspec_kw={"height_ratios": [1, 0.75]}, constrained_layout=True)
    vmin, vmax = min(np.nanmin(base), np.nanmin(opt)), max(np.nanmax(base), np.nanmax(opt))
    first = axes[0, 0].imshow(base, aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax)
    axes[0, 0].set_title("DBP35 experimental effect")
    axes[0, 1].imshow(opt, aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax)
    axes[0, 1].set_title("DBP35opt experimental effect")
    limit = np.nanmax(np.abs(delta))
    third = axes[1, 0].imshow(delta, aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
    axes[1, 0].set_title("DBP35opt - DBP35 effect")
    colors = np.where(position.sensitivity_change >= 0, "#B42318", "#175CD3")
    axes[1, 1].bar(position.position, position.sensitivity_change, color=colors)
    axes[1, 1].axhline(0, color="black", lw=0.8)
    axes[1, 1].set_title("Position sensitivity change")
    axes[1, 1].set_xlabel("Position")
    axes[1, 1].set_ylabel("Delta mean absolute effect")
    for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
        ax.set_xticks(range(14), range(1, 15))
        ax.set_yticks(range(4), order)
        ax.set_xlabel("Position")
        ax.set_ylabel("Mutant base")
    fig.colorbar(first, ax=axes[0, :], shrink=0.82, pad=0.02, label="Normalized effect (-PE/FITC)")
    fig.colorbar(third, ax=axes[1, 0], shrink=0.88, pad=0.03, label="Effect shift")
    fig.suptitle("DBP35 to DBP35opt experimental perturbation (assay contexts differ)")
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"Figure4_DBP35_DBP35opt.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_holdout_report(per: pd.DataFrame, summary: pd.DataFrame, frozen: dict) -> None:
    med = summary.set_index("metric")["median"]
    nulls = pd.read_csv(RESULTS / "holdout_permutation_nulls.tsv", sep="\t")
    lines = [
        "# v1.3.1 Prospective Holdout Evaluation", "",
        "Prediction freeze commit: `ff68394`. Labels were joined only after that commit. The prediction artifact's cached SHA256 is `" + frozen["hash"] + "`; Stage 2 checked only its cached path, size, and mtime and did not re-hash it.", "",
        "## Per-protein decomposition", "",
        "| protein | global rho | position rho | identity residual rho | identity pairwise accuracy | identity permutation p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in per.iterrows():
        p = nulls[(nulls.protein_id == row.protein_id) & (nulls.metric == "within_pairwise_accuracy")].empirical_p.iloc[0]
        lines.append(f"| {row.protein_id} | {row.global_spearman:.3f} | {row.position_spearman:.3f} | {row.within_residual_spearman:.3f} | {row.within_pairwise_accuracy:.3f} | {p:.3f} |")
    lines += [
        "", "## Three-protein summary", "",
        f"- Global median Spearman: {med['global_spearman']:.3f}.",
        f"- Position-sensitivity median Spearman: {med['position_spearman']:.3f}.",
        f"- Within-position identity median Spearman: {med['within_residual_spearman']:.3f}.",
        f"- Within-position pairwise median accuracy: {med['within_pairwise_accuracy']:.3f} (chance 0.5).",
        "- The protein-level sample size is three. Bootstrap intervals are descriptive and no mutation-row pseudo-replication is used for biological inference.",
        "", "## Interpretation gate", "",
    ]
    position_median, identity_median, pair_median = med["position_spearman"], med["within_residual_spearman"], med["within_pairwise_accuracy"]
    if position_median > identity_median and abs(pair_median - 0.5) <= 0.1:
        lines.append("The prospective direction is consistent with Scenario D: position performance exceeds residual identity performance and identity pairwise accuracy remains near chance. This supports the failure-mode decomposition, but n=3 is too small for a broad significance claim.")
    elif identity_median > position_median:
        lines.append("The prospective direction does not reproduce the development position-over-identity pattern (Scenario C-like). The conclusion must be downgraded and protein heterogeneity examined without changing metrics.")
    else:
        lines.append("The prospective result is mixed and does not satisfy a pre-specified model-development gate.")
    lines += [
        "", "Competition replicate identity reproducibility remains unavailable because the official workbook contains published means only. Therefore this holdout result alone cannot distinguish model-specific identity failure from a weakly identifiable assay endpoint, and v1.4 GNN development is not yet justified.",
    ]
    (REPORTS / "V1_3_HOLDOUT_EVALUATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse-evaluation", action="store_true", help="Reuse completed metric/permutation tables and only regenerate reports/figures.")
    args = parser.parse_args()
    FIGURES.mkdir(parents=True, exist_ok=True)
    frozen = check_freeze()
    if args.reuse_evaluation:
        stored = pd.read_csv(RESULTS / "holdout_decomposition.tsv", sep="\t")
        per = stored[stored.row_type == "protein"].copy()
        summary = pd.read_csv(RESULTS / "holdout_decomposition_summary.tsv", sep="\t")
    else:
        per, summary, _ = evaluate_holdout()
    figure2_revised(per)
    figure3_reproducibility(per)
    figure4_dbp35opt()
    write_holdout_report(per, summary, frozen)
    print(json.dumps({
        "status": "complete_after_freeze_commit",
        "freeze_commit": "ff68394",
        "per_protein": per[["protein_id", "global_spearman", "position_spearman", "within_residual_spearman", "within_pairwise_accuracy"]].to_dict("records"),
        "summary": summary.to_dict("records"),
    }, indent=2))


if __name__ == "__main__":
    main()
