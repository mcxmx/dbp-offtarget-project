from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
DOCS = ensure_dir(ROOT / "docs" / "v0_4")
RESULTS = ensure_dir(ROOT / "results" / "v0_4")
TABLES = ensure_dir(RESULTS / "tables")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4")


def fmt(value: float | int | str | None, digits: int = 3) -> str:
    if value is None:
        return "NA"
    try:
        if pd.isna(value):
            return "NA"
    except TypeError:
        pass
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}f}"
    return str(value)


def macro_spearman(macro: pd.DataFrame, baseline: str) -> float:
    row = macro[(macro["baseline"] == baseline) & (macro["metric"] == "spearman")]
    if row.empty:
        return np.nan
    return float(row.iloc[0]["median"])


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    out = df.copy()
    if max_rows is not None:
        out = out.head(max_rows).copy()
    for col in out.columns:
        out[col] = out[col].map(lambda value: fmt(value) if isinstance(value, (float, np.floating)) else ("NA" if pd.isna(value) else str(value)))
    headers = list(out.columns)
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in out.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in headers) + " |")
    return "\n".join(lines)


def write_gap_analysis() -> None:
    macro = pd.read_csv(TABLES / "baseline_performance_macro.csv")
    per = pd.read_csv(TABLES / "baseline_performance_per_protein.csv")
    gap = pd.read_csv(TABLES / "baseline_gap_summary.csv")
    resolution = pd.read_csv(TABLES / "failure_resolution_summary.csv")
    failures = pd.read_parquet(TABLES / "baseline_failure_cases.parquet")
    noise = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "experimental_noise_ceiling.csv")
    distance = pd.read_csv(TABLES / "performance_by_sequence_distance.csv")

    replicate_ref = float(noise.loc[noise["score_type"] == "e_score", "spearman_correlation"].median())
    best_sequence = (
        gap[gap["baseline"].str.startswith("sequence_")]
        .sort_values("macro_median_spearman", ascending=False)
        .iloc[0]
    )
    nampnn = gap[gap["baseline"] == "NA-MPNN_structural_ppm"].iloc[0]
    total_disagreement = int(gap["disagreement_candidates_total"].dropna().iloc[0])
    resolved = resolution[resolution["baseline"] == "NA-MPNN_structural_ppm"]["n_resolved"].sum()
    evaluable = resolution[resolution["baseline"] == "NA-MPNN_structural_ppm"]["n_evaluable_candidates"].sum()
    high_low = int((failures["failure_category"] == "experimental_high_sequence_proxy_low").sum())
    low_high = int((failures["failure_category"] == "sequence_proxy_high_experimental_low").sum())

    hardest_sequence = (
        per[(per["baseline"] == best_sequence["baseline"]) & per["spearman"].notna()]
        .sort_values("spearman")
        .iloc[0]
    )
    seq_distance = distance[(distance["baseline"] == best_sequence["baseline"]) & distance["spearman"].notna()]
    hardest_distance = seq_distance.sort_values("spearman").iloc[0] if not seq_distance.empty else None

    per_table = markdown_table(per[["protein_id", "baseline", "spearman", "evaluation_status"]])
    resolution_table = markdown_table(resolution)
    gap_table = markdown_table(gap)

    text = f"""# v0.4 Baseline Gap Analysis

Analysis date: 2026-09-02

## Summary

The fixed benchmark is v0.3.1 GSE237017 designed-DBP uPBM: 7 proteins and 57,344 protein-RC-class experimental units. Scores are processed experimental PBM E-score consensus values and are evaluated as per-protein DNA ranking targets.

Best sequence-only baseline: `{best_sequence['baseline']}` with macro median Spearman {fmt(best_sequence['macro_median_spearman'])}. The empirical replicate agreement reference is median E-score replicate Spearman {fmt(replicate_ref)}. The gap from the best sequence-only median to this reference is {fmt(best_sequence['gap_to_reference'])}.

DeepPBS was not fairly evaluable in the current Windows environment because the official preprocessing stack requires additional structure-processing dependencies. SimpleProteinConditionalBaseline is intentionally untrained because no assay-matched natural PBM/uPBM training set has been added yet.

NA-MPNN ran only as a structural diagnostic for DBP35 and DBP48. Its macro median Spearman over these two evaluable proteins was {fmt(nampnn['macro_median_spearman'])}. DBP48/8TAC appears in NA-MPNN split files and is not a zero-shot result.

## What Current Baselines Can Explain

- Sequence-only 3-mer similarity is the strongest Tier 0 baseline, but remains well below empirical replicate agreement.
- DBP35 NA-MPNN diagnostic Spearman is positive, but DBP48 is negative despite having an experimental structure and a training/validation overlap warning.
- Coverage is the main structure-aware bottleneck: five of seven designed DBPs lack a public structure/model found in the checked official sources.

## Disagreement Cases

v0.3.1 defined {total_disagreement} sequence-vs-experiment disagreement candidates. In v0.4, NA-MPNN predictions exist for {int(evaluable)} of these candidates and rank {int(resolved)} of them in the top 10% of NA-MPNN scores. This is a diagnostic resolution rate of {fmt(resolved / evaluable if evaluable else np.nan)} among NA-MPNN-evaluable disagreement candidates, not across the whole benchmark.

The v0.4 failure table contains {len(failures)} sequence-vs-experiment ranking cases: {high_low} experimental-high/sequence-low cases and {low_high} sequence-high/experimental-low cases.

## Hardest Observed Regimes

For the best sequence-only baseline (`{best_sequence['baseline']}`), the lowest per-protein Spearman is {fmt(hardest_sequence['spearman'])} on {hardest_sequence['protein_id']}. The lowest evaluated motif-distance stratum is {hardest_distance['protein_id']} distance {hardest_distance['motif_distance_bin']} with Spearman {fmt(hardest_distance['spearman'])} if using the same baseline.

## Per-Protein Metrics

{per_table}

## Macro Gap Table

{gap_table}

## Disagreement Resolution Table

{resolution_table}

## Interpretation Limits

These results do not show that a new model is better than DeepPBS or NA-MPNN. They show that the current benchmark exposes a large sequence-only gap and that existing structure-aware methods are not yet comprehensively evaluable on the seven designed DBPs with public structures available in this repository.
"""
    (RESULTS / "BASELINE_GAP_ANALYSIS.md").write_text(text, encoding="utf-8")


def write_go_no_go() -> None:
    gap = pd.read_csv(TABLES / "baseline_gap_summary.csv")
    per = pd.read_csv(TABLES / "baseline_performance_per_protein.csv")
    overlap = pd.read_csv(METADATA / "baseline_data_overlap_audit.csv")
    structures = pd.read_csv(METADATA / "designed_dbp_structure_manifest.csv")
    noise = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "experimental_noise_ceiling.csv")

    rep = float(noise.loc[noise["score_type"] == "e_score", "spearman_correlation"].median())
    best_sequence = (
        gap[gap["baseline"].str.startswith("sequence_")]
        .sort_values("macro_median_spearman", ascending=False)
        .iloc[0]
    )
    nampnn_n = int(gap.loc[gap["baseline"] == "NA-MPNN_structural_ppm", "n_proteins_evaluated"].iloc[0])
    not_found = int((structures["structure_confidence"].fillna("") == "none").sum())
    high_overlap = overlap[overlap["risk_level"].astype(str).str.contains("high", case=False, na=False)]

    text = f"""# v0.4 New Model Go / No-Go

Decision: CONDITIONAL GO

## Rationale

The fixed v0.3.1 benchmark is usable for baseline arena work, and the sequence-only gap is clear: best sequence-only macro median Spearman is {fmt(best_sequence['macro_median_spearman'])}, versus empirical replicate agreement reference {fmt(rep)}.

However, this is not a STRONG GO yet. DeepPBS was not fairly runnable in the current environment, and NA-MPNN was evaluable for only {nampnn_n}/7 designed DBPs. Five designed DBPs have no public structure/model found in the checked official sources. DBP48/8TAC also has a high overlap risk because `8tac` appears in NA-MPNN split files.

## Required Before Strong Claims

- Run DeepPBS in a supported Linux/container environment or explicitly document it as not comparable.
- Expand structure availability or use structure-free baselines so all seven designed DBPs can be evaluated.
- Add assay-matched natural PBM/uPBM training data before training the Tier 1 protein-conditioned baseline.
- Keep DBP48/8TAC separate from zero-shot claims because of the detected NA-MPNN split overlap.

## Gate Interpretation

CONDITIONAL GO means v0.4 should continue into stronger baseline work and a careful Tier 1 assay-matched training setup. It does not authorize claims that existing strong baselines have systematically failed across the full designed-DBP benchmark.

## Key Constraints

- No new final model should be proposed until baseline coverage is improved.
- Any future comparison must remain per-protein and RC-class grouped.
- uPBM E-score remains an experimental specificity ranking signal, not affinity or in vivo binding.
"""
    (RESULTS / "NEW_MODEL_GO_NO_GO.md").write_text(text, encoding="utf-8")


def write_model_requirements() -> None:
    gap = pd.read_csv(TABLES / "baseline_gap_summary.csv")
    resolution = pd.read_csv(TABLES / "failure_resolution_summary.csv")
    best_sequence = (
        gap[gap["baseline"].str.startswith("sequence_")]
        .sort_values("macro_median_spearman", ascending=False)
        .iloc[0]
    )
    resolved = resolution[resolution["baseline"] == "NA-MPNN_structural_ppm"]["n_resolved"].sum()
    evaluable = resolution[resolution["baseline"] == "NA-MPNN_structural_ppm"]["n_evaluable_candidates"].sum()
    text = f"""# v0.4 Proposed Model Requirements

This document is requirements-only. v0.4 does not implement a final model.

## Requirements Driven by Current Failure Analysis

1. The future baseline arena needs assay-matched natural PBM/uPBM data before training a protein-conditioned model. Without that control, natural-to-designed performance drops would be confounded by assay shift.
2. The model must be evaluated as a per-protein RC-class ranking task. The best current sequence-only baseline is `{best_sequence['baseline']}` with macro median Spearman {fmt(best_sequence['macro_median_spearman'])}, below empirical replicate agreement.
3. The method should not require a public complex structure for every designed DBP unless structure generation/mapping is itself benchmarked. v0.4 found public structure support for only DBP35 and DBP48 in checked sources.
4. The model must report missing predictions as not evaluable. Missing DeepPBS or SimpleProteinConditionalBaseline scores must never be filled with zero.
5. The model should explicitly target sequence-vs-experiment disagreement candidates. NA-MPNN resolved {int(resolved)}/{int(evaluable)} evaluable v0.3.1 disagreement candidates in this diagnostic setup.
6. DBP48 must be handled as a non-zero-shot diagnostic for NA-MPNN because 8TAC appears in NA-MPNN split files.

## Not Yet Justified

- A large target-anchored neural architecture.
- Cross-protein absolute E-score calibration.
- Claims that designed DBPs are OOD without assay-matched natural PBM/uPBM controls.
"""
    (DOCS / "PROPOSED_MODEL_REQUIREMENTS.md").write_text(text, encoding="utf-8")


def write_progress() -> None:
    macro = pd.read_csv(TABLES / "baseline_performance_macro.csv")
    per = pd.read_csv(TABLES / "baseline_performance_per_protein.csv")
    resolution = pd.read_csv(TABLES / "failure_resolution_summary.csv")
    structures = pd.read_csv(METADATA / "designed_dbp_structure_manifest.csv")
    overlap = pd.read_csv(METADATA / "baseline_data_overlap_audit.csv")
    scored = pd.read_parquet(ROOT / "data" / "processed" / "v0_4" / "v0_4_scored_candidates.parquet")

    best_seq = macro[(macro["metric"] == "spearman") & macro["baseline"].str.startswith("sequence_")].sort_values("median", ascending=False).iloc[0]
    nampnn_med = macro_spearman(macro, "NA-MPNN_structural_ppm")
    evaluated_nampnn = per[(per["baseline"] == "NA-MPNN_structural_ppm") & (per["evaluation_status"] == "evaluated")]["protein_id"].tolist()
    high_overlap = overlap[overlap["risk_level"].astype(str).str.contains("high", case=False, na=False)]
    resolved = int(resolution[resolution["baseline"] == "NA-MPNN_structural_ppm"]["n_resolved"].sum())
    evaluable = int(resolution[resolution["baseline"] == "NA-MPNN_structural_ppm"]["n_evaluable_candidates"].sum())
    n_rc = scored[["protein_id", "canonical_7mer"]].drop_duplicates().shape[0]
    no_structure = structures.loc[structures["structure_confidence"].fillna("") == "none", "protein_id"].tolist()

    text = f"""# v0.4 进展报告

本轮在 v0.3.1 已通过 validation 的 GSE237017 designed-DBP uPBM benchmark 上，建立了 strong baseline arena 的第一版。v0.3.1 数据保持冻结，v0.4 新结果单独保存在 `data/processed/v0_4/`、`metadata/v0_4/` 和 `results/v0_4/`。

已完成 DeepPBS、NA-MPNN 和 Tier 1 SimpleProteinConditionalBaseline 的可行性审计。DeepPBS 本轮未公平运行，原因是官方预处理依赖 Linux/结构特征工具链和额外图神经网络依赖；SimpleProteinConditionalBaseline 只保留为未训练的 protein-conditioned 接口，因为尚未加入 assay-matched natural PBM/uPBM 训练集。NA-MPNN 使用官方 specificity checkpoint 成功对 DBP35 和 DBP48 产生结构 PPM 诊断预测，但 DBP48/8TAC 在 NA-MPNN split 文件中出现，不能作为 zero-shot 结果。

当前评估覆盖 {n_rc} 个 protein-RC-class 单位。最好的 sequence-only baseline 是 `{best_seq['baseline']}`，macro median Spearman 为 {fmt(best_seq['median'])}。NA-MPNN 诊断结果覆盖 {len(evaluated_nampnn)}/7 个 DBP，macro median Spearman 为 {fmt(nampnn_med)}。五个暂未能结构评估的 DBP 是：{', '.join(no_structure)}。v0.3.1 的 1,515 个 sequence-vs-experiment disagreement candidates 中，NA-MPNN 可评估 {evaluable} 个，按预设 top-10% 规则解析 {resolved} 个。

本轮生成了 baseline performance、per-protein heatmap、prediction-vs-experiment、replicate reference、motif-distance 分层和 failure landscape 六张图。最终 gate 为 `CONDITIONAL GO`：sequence-only gap 明显，但 DeepPBS/NA-MPNN 的全覆盖强 baseline 证据仍不足。下一步最应优先补充 assay-matched natural PBM/uPBM training control，并在可复现环境中补跑 DeepPBS 或替代结构-aware baseline。
"""
    (RESULTS / "V0_4_PROGRESS.md").write_text(text, encoding="utf-8")


def main() -> None:
    write_gap_analysis()
    write_go_no_go()
    write_model_requirements()
    write_progress()
    print("Wrote v0.4 reports")


if __name__ == "__main__":
    main()
