from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.structure_baselines import score_structural_ppm
from src.utils import ensure_dir, project_root
from src.v0_4_evaluation import bootstrap_metric_ci, compute_ranking_metrics, macro_summary


ROOT = project_root()
DOCS = ensure_dir(ROOT / "docs" / "v0_4")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4")
DATA = ensure_dir(ROOT / "data" / "processed" / "v0_4")
RESULTS = ensure_dir(ROOT / "results" / "v0_4")
TABLES = ensure_dir(RESULTS / "tables")
FIGURES = ensure_dir(RESULTS / "figures")
EXTERNAL = ensure_dir(RESULTS / "external_runs")


def load_rc_benchmark() -> pd.DataFrame:
    df = pd.read_parquet(ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet").copy()
    df["experimental_percentile_within_protein"] = df.groupby("protein_id")["experimental_escore_consensus"].rank(pct=True, ascending=True)
    return df


def load_sequence_baseline() -> pd.DataFrame:
    df = pd.read_parquet(ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1.parquet").copy()
    df = df.rename(
        columns={
            "hamming_similarity_to_paper_motif_rc_aware": "sequence_hamming_score",
            "edit_similarity_to_paper_motif_rc_aware": "sequence_edit_score",
            "kmer3_jaccard_to_paper_motif_rc_aware": "sequence_kmer3_score",
            "kmer4_jaccard_to_paper_motif_rc_aware": "sequence_kmer4_score",
        }
    )
    return df


def score_external_structural_baselines() -> pd.DataFrame:
    rows = []
    # DBP35 theoretical structure from official dbp_design repo.
    dbp35_npz = ROOT / "results" / "v0_4" / "external_runs" / "nampnn_dbp035" / "specificity" / "DBP035.npz"
    if dbp35_npz.exists():
        rows.append(
            score_structural_ppm(
                dbp35_npz,
                protein_id="DBP35",
                structure_id="dbp_design:DBP035.pdb",
                model_version="NA-MPNN specificity s_70114.pt",
                prediction_type="partial_structural_ppm_best_window_log_probability",
            )
        )
    dbp48_npz = ROOT / "results" / "v0_4" / "external_runs" / "nampnn_8tac" / "specificity" / "8TAC.npz"
    if dbp48_npz.exists():
        rows.append(
            score_structural_ppm(
                dbp48_npz,
                protein_id="DBP48",
                structure_id="RCSB:8TAC",
                model_version="NA-MPNN specificity s_70114.pt",
                prediction_type="partial_structural_ppm_best_window_log_probability",
            )
        )
    if not rows:
        return pd.DataFrame(columns=[
            "protein_id",
            "canonical_7mer",
            "oriented_7mer",
            "reverse_complement_7mer",
            "prediction_score",
            "prediction_type",
            "structure_id",
            "model_version",
            "best_chain_label",
            "best_window_start",
            "best_window_end",
            "best_window_sequence",
            "prediction_orientation",
            "source_npz",
        ])
    return pd.concat(rows, ignore_index=True)


def add_prediction_columns(rc: pd.DataFrame, seq: pd.DataFrame, structural: pd.DataFrame) -> pd.DataFrame:
    df = rc.merge(
        seq[[
            "protein_id",
            "canonical_7mer",
            "sequence_hamming_score",
            "sequence_edit_score",
            "sequence_kmer3_score",
            "sequence_kmer4_score",
        ]],
        on=["protein_id", "canonical_7mer"],
        how="left",
    )
    if not structural.empty:
        df = df.merge(
            structural[["protein_id", "canonical_7mer", "prediction_score"]].rename(
                columns={"prediction_score": "structural_ppm_score"}
            ),
            on=["protein_id", "canonical_7mer"],
            how="left",
        )
    else:
        df["structural_ppm_score"] = np.nan
    return df


def per_protein_metrics(scored: pd.DataFrame) -> pd.DataFrame:
    baselines = {
        "sequence_hamming": "sequence_hamming_score",
        "sequence_edit": "sequence_edit_score",
        "sequence_kmer3": "sequence_kmer3_score",
        "sequence_kmer4": "sequence_kmer4_score",
        "NA-MPNN_structural_ppm": "structural_ppm_score",
        "DeepPBS": "deeppbs_score",
        "SimpleProteinConditionalBaseline": "simple_pc_score",
    }
    rows = []
    for protein_id, group in scored.groupby("protein_id"):
        for baseline, score_col in baselines.items():
            sub = group[["experimental_escore_consensus", score_col, "canonical_7mer"]].dropna()
            if sub.empty or sub[score_col].notna().sum() == 0:
                if baseline == "NA-MPNN_structural_ppm":
                    status = "not_evaluable_missing_public_structure"
                elif baseline == "DeepPBS":
                    status = "not_evaluable_current_environment"
                elif baseline == "SimpleProteinConditionalBaseline":
                    status = "not_trained_no_assay_matched_training_data"
                else:
                    status = "no_scores"
                rows.append(
                    {
                        "protein_id": protein_id,
                        "baseline": baseline,
                        "spearman": np.nan,
                        "ndcg_1pct": np.nan,
                        "ndcg_5pct": np.nan,
                        "pairwise_accuracy": np.nan,
                        "top1pct_recovery": np.nan,
                        "n_rc_classes": 0,
                        "evaluation_status": status,
                    }
                )
                continue
            metrics = compute_ranking_metrics(sub.rename(columns={score_col: "prediction_score"}), "experimental_escore_consensus", "prediction_score")
            rows.append(
                {
                    "protein_id": protein_id,
                    "baseline": baseline,
                    "spearman": metrics.spearman,
                    "ndcg_1pct": metrics.ndcg_1pct,
                    "ndcg_5pct": metrics.ndcg_5pct,
                    "pairwise_accuracy": metrics.pairwise_accuracy,
                    "top1pct_recovery": metrics.top1pct_recovery,
                    "n_rc_classes": metrics.n_rc_classes,
                    "evaluation_status": "evaluated",
                }
            )
    return pd.DataFrame(rows)


def protein_bootstrap_ci(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    baselines = {
        "sequence_hamming": "sequence_hamming_score",
        "sequence_edit": "sequence_edit_score",
        "sequence_kmer3": "sequence_kmer3_score",
        "sequence_kmer4": "sequence_kmer4_score",
        "NA-MPNN_structural_ppm": "structural_ppm_score",
        "DeepPBS": "deeppbs_score",
        "SimpleProteinConditionalBaseline": "simple_pc_score",
    }
    for protein_id, group in scored.groupby("protein_id"):
        for baseline, score_col in baselines.items():
            sub = group[["experimental_escore_consensus", score_col, "canonical_7mer"]].dropna()
            if sub.empty or sub[score_col].notna().sum() == 0:
                if baseline == "NA-MPNN_structural_ppm":
                    status = "not_evaluable_missing_public_structure"
                elif baseline == "DeepPBS":
                    status = "not_evaluable_current_environment"
                elif baseline == "SimpleProteinConditionalBaseline":
                    status = "not_trained_no_assay_matched_training_data"
                else:
                    status = "no_scores"
                rows.append(
                    {
                        "protein_id": protein_id,
                        "baseline": baseline,
                        "metric": "spearman",
                        "bootstrap_mean": np.nan,
                        "ci_lower": np.nan,
                        "ci_upper": np.nan,
                        "n_resamples": 0,
                        "evaluation_status": status,
                    }
                )
                continue
            mean_, lo, hi = bootstrap_metric_ci(sub.rename(columns={score_col: "prediction_score"}), "experimental_escore_consensus", "prediction_score", metric="spearman", n_bootstrap=200, seed=42)
            rows.append(
                {
                    "protein_id": protein_id,
                    "baseline": baseline,
                    "metric": "spearman",
                    "bootstrap_mean": mean_,
                    "ci_lower": lo,
                    "ci_upper": hi,
                    "n_resamples": 200,
                    "evaluation_status": "evaluated",
                }
            )
    return pd.DataFrame(rows)


def failure_cases(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein_id, group in scored.groupby("protein_id"):
        g = group.copy()
        g["experimental_percentile"] = g["experimental_escore_consensus"].rank(pct=True, ascending=True)
        g["sequence_percentile"] = g["sequence_kmer3_score"].rank(pct=True, ascending=True)
        for _, row in g.iterrows():
            if row["experimental_percentile"] >= 0.95 and row["sequence_percentile"] <= 0.50:
                category = "experimental_high_sequence_proxy_low"
            elif row["experimental_percentile"] <= 0.50 and row["sequence_percentile"] >= 0.95:
                category = "sequence_proxy_high_experimental_low"
            else:
                continue
            rows.append(
                {
                    "protein_id": protein_id,
                    "canonical_7mer": row["canonical_7mer"],
                    "experimental_escore": row["experimental_escore_consensus"],
                    "experimental_percentile": float(row["experimental_percentile"]),
                    "sequence_score": float(row["sequence_kmer3_score"]),
                    "deeppbs_score": np.nan,
                    "nampnn_score": row["structural_ppm_score"] if pd.notna(row.get("structural_ppm_score")) else np.nan,
                    "simple_pc_score": np.nan,
                    "failure_category": category,
                }
            )
    return pd.DataFrame(rows)


def failure_resolution_summary(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    baseline_cols = {
        "NA-MPNN_structural_ppm": "structural_ppm_score",
        "DeepPBS": "deeppbs_score",
        "SimpleProteinConditionalBaseline": "simple_pc_score",
    }
    for protein_id, group in scored.groupby("protein_id"):
        g = group.copy()
        experimental_threshold = float(g["experimental_escore_consensus"].quantile(0.95))
        sequence_threshold = float(g["sequence_hamming_score"].quantile(0.50))
        candidates = g[
            (g["experimental_escore_consensus"] >= experimental_threshold)
            & (g["sequence_hamming_score"] <= sequence_threshold)
        ].copy()
        for baseline, col in baseline_cols.items():
            evaluable = candidates.dropna(subset=[col]).copy()
            if evaluable.empty:
                rows.append(
                    {
                        "protein_id": protein_id,
                        "baseline": baseline,
                        "n_v0_3_1_disagreement_candidates": int(len(candidates)),
                        "n_evaluable_candidates": 0,
                        "n_resolved": 0,
                        "resolution_rate": np.nan,
                        "resolution_definition": "candidate predicted in top 10% within protein by protein-conditioned baseline score",
                        "evaluation_status": "not_evaluable",
                    }
                )
                continue
            g["prediction_percentile"] = g[col].rank(pct=True, ascending=True)
            evaluable = g.loc[evaluable.index]
            resolved = evaluable[evaluable["prediction_percentile"] >= 0.90]
            rows.append(
                {
                    "protein_id": protein_id,
                    "baseline": baseline,
                    "n_v0_3_1_disagreement_candidates": int(len(candidates)),
                    "n_evaluable_candidates": int(len(evaluable)),
                    "n_resolved": int(len(resolved)),
                    "resolution_rate": float(len(resolved) / len(evaluable)),
                    "resolution_definition": "candidate predicted in top 10% within protein by protein-conditioned baseline score",
                    "evaluation_status": "evaluated",
                }
            )
    return pd.DataFrame(rows)


def performance_by_distance(scored: pd.DataFrame) -> pd.DataFrame:
    df = scored.copy()
    df["motif_length"] = df["designed_binding_site_motif"].str.len()
    df["motif_hamming_distance"] = ((1.0 - df["sequence_hamming_score"]) * df["motif_length"]).round().astype(int)
    df["motif_distance_bin"] = df["motif_hamming_distance"].map(lambda x: str(x) if x < 3 else "3+")
    baselines = {
        "sequence_hamming": "sequence_hamming_score",
        "sequence_edit": "sequence_edit_score",
        "sequence_kmer3": "sequence_kmer3_score",
        "sequence_kmer4": "sequence_kmer4_score",
        "NA-MPNN_structural_ppm": "structural_ppm_score",
        "DeepPBS": "deeppbs_score",
        "SimpleProteinConditionalBaseline": "simple_pc_score",
    }
    rows = []
    for (protein_id, distance_bin), group in df.groupby(["protein_id", "motif_distance_bin"]):
        for baseline, col in baselines.items():
            sub = group[["experimental_escore_consensus", col, "canonical_7mer"]].dropna()
            if len(sub) < 20:
                rows.append(
                    {
                        "protein_id": protein_id,
                        "motif_distance_bin": distance_bin,
                        "baseline": baseline,
                        "spearman": np.nan,
                        "ndcg_5pct": np.nan,
                        "n_rc_classes": int(len(sub)),
                        "evaluation_status": "too_few_or_no_scores",
                    }
                )
                continue
            metrics = compute_ranking_metrics(sub.rename(columns={col: "prediction_score"}), "experimental_escore_consensus", "prediction_score")
            rows.append(
                {
                    "protein_id": protein_id,
                    "motif_distance_bin": distance_bin,
                    "baseline": baseline,
                    "spearman": metrics.spearman,
                    "ndcg_5pct": metrics.ndcg_5pct,
                    "n_rc_classes": metrics.n_rc_classes,
                    "evaluation_status": "evaluated",
                }
            )
    return pd.DataFrame(rows)


def performance_by_group(per_protein: pd.DataFrame) -> pd.DataFrame:
    groups = pd.read_csv(ROOT / "metadata" / "v0_3_1" / "designed_dbp_target_groups.csv")
    merged = per_protein.merge(groups, on="protein_id", how="left")
    rows = []
    for group_col in ["protein_sequence_cluster", "original_target_group", "assay_target_group", "motif_group"]:
        for (baseline, group_id), group in merged.groupby(["baseline", group_col]):
            vals = group["spearman"].dropna().to_numpy(dtype=float)
            rows.append(
                {
                    "group_type": group_col,
                    "group_id": group_id,
                    "baseline": baseline,
                    "n_proteins_with_metric": int(vals.size),
                    "median_spearman": float(np.median(vals)) if vals.size else np.nan,
                    "mean_spearman": float(np.mean(vals)) if vals.size else np.nan,
                }
            )
    return pd.DataFrame(rows)


def build_gap_summary(per_protein: pd.DataFrame, noise_path: Path, disagreement_count: int) -> pd.DataFrame:
    noise = pd.read_csv(noise_path)
    rep = noise[noise["score_type"] == "e_score"]["spearman_correlation"].median()
    rows = []
    for baseline, group in per_protein.groupby("baseline"):
        vals = group["spearman"].dropna()
        rows.append(
            {
                "baseline": baseline,
                "macro_median_spearman": float(vals.median()) if not vals.empty else np.nan,
                "macro_mean_spearman": float(vals.mean()) if not vals.empty else np.nan,
                "n_proteins_evaluated": int(vals.size),
                "replicate_reference_spearman": float(rep),
                "gap_to_reference": float(rep - vals.median()) if not vals.empty else np.nan,
                "disagreement_candidates_total": disagreement_count,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    rc = load_rc_benchmark()
    seq = load_sequence_baseline()
    structural = score_external_structural_baselines()
    if not structural.empty:
        structural["source_npz"] = structural["source_npz"].map(
            lambda value: str(Path(value).resolve().relative_to(ROOT)) if Path(value).resolve().is_relative_to(ROOT) else value
        )
    structural.to_parquet(DATA / "nampnn_structural_predictions.parquet", index=False)
    structural.to_csv(TABLES / "nampnn_predictions.csv", index=False)
    structural.to_parquet(DATA / "nampnn_predictions.parquet", index=False)
    structural.to_parquet(TABLES / "nampnn_predictions.parquet", index=False)
    empty_prediction = pd.DataFrame(
        columns=["protein_id", "canonical_7mer", "prediction_score", "prediction_type", "structure_id", "model_version"]
    )
    empty_prediction.to_parquet(DATA / "deeppbs_predictions.parquet", index=False)
    empty_prediction.to_parquet(TABLES / "deeppbs_predictions.parquet", index=False)
    empty_prediction.to_parquet(DATA / "simple_pc_predictions.parquet", index=False)
    scored = add_prediction_columns(rc, seq, structural)
    scored["deeppbs_score"] = np.nan
    scored["simple_pc_score"] = np.nan
    scored.to_parquet(DATA / "v0_4_scored_candidates.parquet", index=False)
    scored.head(500).to_csv(TABLES / "v0_4_scored_candidates_preview.csv", index=False)

    per_protein = per_protein_metrics(scored)
    per_protein.to_csv(TABLES / "baseline_performance_per_protein.csv", index=False)

    macro_rows = []
    for baseline, group in per_protein.groupby("baseline"):
        for metric in ["spearman", "ndcg_1pct", "ndcg_5pct", "pairwise_accuracy", "top1pct_recovery"]:
            vals = group[metric].dropna().to_numpy(dtype=float)
            macro_rows.append(
                {
                    "baseline": baseline,
                    "metric": metric,
                    "mean": float(np.mean(vals)) if vals.size else np.nan,
                    "median": float(np.median(vals)) if vals.size else np.nan,
                    "std": float(np.std(vals, ddof=1)) if vals.size > 1 else np.nan,
                    "n_proteins": int(vals.size),
                }
            )
    pd.DataFrame(macro_rows).to_csv(TABLES / "baseline_performance_macro.csv", index=False)
    bootstrap = protein_bootstrap_ci(scored)
    macro_bootstrap_rows = []
    rng = np.random.default_rng(42)
    for baseline, group in per_protein.groupby("baseline"):
        vals = group["spearman"].dropna().to_numpy(dtype=float)
        if vals.size:
            boot = [float(np.median(rng.choice(vals, size=vals.size, replace=True))) for _ in range(200)]
            macro_bootstrap_rows.append(
                {
                    "protein_id": "macro_protein_bootstrap",
                    "baseline": baseline,
                    "metric": "spearman",
                    "bootstrap_mean": float(np.mean(boot)),
                    "ci_lower": float(np.percentile(boot, 2.5)),
                    "ci_upper": float(np.percentile(boot, 97.5)),
                    "n_resamples": 200,
                    "evaluation_status": "evaluated",
                }
            )
        else:
            macro_bootstrap_rows.append(
                {
                    "protein_id": "macro_protein_bootstrap",
                    "baseline": baseline,
                    "metric": "spearman",
                    "bootstrap_mean": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "n_resamples": 0,
                    "evaluation_status": "not_evaluable",
                }
            )
    bootstrap = pd.concat([bootstrap, pd.DataFrame(macro_bootstrap_rows)], ignore_index=True)
    bootstrap.to_csv(TABLES / "baseline_bootstrap_ci.csv", index=False)

    failures = failure_cases(scored)
    failures.to_parquet(TABLES / "baseline_failure_cases.parquet", index=False)
    failures.to_csv(TABLES / "baseline_failure_cases.csv", index=False)
    failure_resolution_summary(scored).to_csv(TABLES / "failure_resolution_summary.csv", index=False)
    performance_by_distance(scored).to_csv(TABLES / "performance_by_sequence_distance.csv", index=False)
    performance_by_group(per_protein).to_csv(TABLES / "baseline_performance_by_group.csv", index=False)

    disagreement = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "all_disagreement_candidate_counts.csv")
    gap = build_gap_summary(per_protein, ROOT / "results" / "v0_3_1" / "tables" / "experimental_noise_ceiling.csv", int(disagreement["n_disagreement"].sum()))
    gap.to_csv(TABLES / "baseline_gap_summary.csv", index=False)

    # v0.4.1-style provenance files required by tests and reports.
    per_protein.to_csv(TABLES / "baseline_performance_per_protein_full.csv", index=False)
    print(per_protein.to_string(index=False))
    print(gap.to_string(index=False))


if __name__ == "__main__":
    main()
