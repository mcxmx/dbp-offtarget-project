from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
METADATA_DIR = ensure_dir(ROOT / "metadata")
RESULTS_DIR = ensure_dir(ROOT / "results")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")


def count_value(df: pd.DataFrame, column: str, value: str) -> int:
    return int((df[column] == value).sum())


def quality_summary(curation: pd.DataFrame, pairs_v02: pd.DataFrame, benchmark: pd.DataFrame, experimental: pd.DataFrame) -> dict[str, int]:
    return {
        "raw_pdb_pairs": int(len(curation)),
        "curated_core_pairs": int(len(pairs_v02)),
        "sequence_specific_pairs": count_value(curation, "sequence_specificity_class", "sequence_specific"),
        "weakly_sequence_specific_pairs": count_value(curation, "sequence_specificity_class", "weakly_sequence_specific"),
        "non_specific_pairs": count_value(curation, "sequence_specificity_class", "non_specific"),
        "guide_dependent_pairs": count_value(curation, "sequence_specificity_class", "guide_dependent"),
        "lesion_specific_pairs": count_value(curation, "sequence_specificity_class", "lesion_specific"),
        "structure_specific_pairs": count_value(curation, "sequence_specificity_class", "structure_specific"),
        "designed_sequence_specific_pairs": count_value(curation, "sequence_specificity_class", "designed_sequence_specific"),
        "uncertain_pairs": count_value(curation, "sequence_specificity_class", "uncertain"),
        "pdb_quantitative_ground_truth_pairs": int(curation["has_quantitative_specificity_ground_truth"].astype(str).str.lower().eq("true").sum()),
        "pdb_without_quantitative_ground_truth_pairs": int(curation["has_quantitative_specificity_ground_truth"].astype(str).str.lower().ne("true").sum()),
        "benchmark_v0_2_rows": int(len(benchmark)),
        "single_mutants_v0_2": int((benchmark["candidate_type"] == "single_mut").sum()),
        "double_mutants_v0_2": int((benchmark["candidate_type"] == "double_mut").sum()),
        "random_negatives_v0_2": int(benchmark["candidate_type"].isin(["gc_matched_random", "random_dna"]).sum()),
        "experimental_specificity_proteins": int(experimental["protein_id"].nunique()),
        "experimental_specificity_rows": int(len(experimental)),
    }


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "(no rows)"
    columns = [str(col) for col in df.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in df.columns) + " |")
    return "\n".join(lines)


def write_quality_report(summary: dict[str, int], curation: pd.DataFrame, benchmark_counts: pd.DataFrame) -> None:
    class_counts = (
        curation["sequence_specificity_class"]
        .value_counts()
        .rename_axis("sequence_specificity_class")
        .reset_index(name="count")
    )
    use_counts = (
        curation["recommended_use"]
        .value_counts()
        .rename_axis("recommended_use")
        .reset_index(name="count")
    )
    bench_md = markdown_table(benchmark_counts)
    report = f"""# Benchmark Quality Report v0.2

Generated from repository tables on 2026-08-29.

## Summary Counts

- Raw historical PDB pairs: {summary['raw_pdb_pairs']}
- Curated v0.2 core benchmark pairs: {summary['curated_core_pairs']}
- Sequence-specific pairs: {summary['sequence_specific_pairs']}
- Non-specific pairs: {summary['non_specific_pairs']}
- Guide-dependent pairs: {summary['guide_dependent_pairs']}
- Lesion-specific pairs: {summary['lesion_specific_pairs']}
- Designed sequence-specific pairs: {summary['designed_sequence_specific_pairs']}
- Uncertain pairs: {summary['uncertain_pairs']}
- PDB pairs with quantitative specificity ground truth: {summary['pdb_quantitative_ground_truth_pairs']}
- PDB pairs without quantitative specificity ground truth: {summary['pdb_without_quantitative_ground_truth_pairs']}
- Experimental specificity pilot proteins: {summary['experimental_specificity_proteins']}
- Experimental specificity pilot rows: {summary['experimental_specificity_rows']}

## Sequence Specificity Classes

{markdown_table(class_counts)}

## Recommended Use

{markdown_table(use_counts)}

## Benchmark v0.2 Candidate Counts

{bench_md}

## v0.1 to v0.2 Corrections

- Split PDB structural evidence from quantitative specificity ground truth.
- Set all PDB-only records to `has_quantitative_specificity_ground_truth=False`.
- Replaced the longest-chain benchmark assumption with curated chain annotation
  plus chain-contact evidence where available.
- Moved guide-dependent, lesion-specific, non-specific, and transposase/substrate
  cases out of the v0.2 core specificity benchmark.
- Replotted v0.2 figures with sequence-only proxy terminology.
- Added an explicit positional-bias analysis for k-mer and combined proxy
  metrics.

## Current Interpretation

The v0.2 structural/mutation benchmark is suitable for reproducible pipeline
testing and sequence-only sanity checks. It is not yet a calibrated
protein-conditioned off-target predictor. The experimental specificity pilot is
a separate Layer C resource based on JASPAR PFM-derived PWM scores, not raw PBM
or HT-SELEX enrichment.
"""
    (RESULTS_DIR / "BENCHMARK_QUALITY_REPORT.md").write_text(report, encoding="utf-8")


def write_weekly_progress(summary: dict[str, int]) -> None:
    report = f"""# Weekly Progress v0.2

本周对 DBP off-target benchmark 原型完成了一次 scientific audit 和 v0.2 修正。首先检查了 v0.1 的数据生成逻辑，发现 PDB complex 曾被误标为 specificity ground truth，并且使用“最长蛋白链/最长 DNA 链”作为 benchmark 选择规则。v0.2 将证据拆分为 structural cognate、direct DNA-binding evidence、sequence-specificity evidence 和 quantitative specificity ground truth；当前 {summary['raw_pdb_pairs']} 个 PDB pair 均不再被视为 quantitative ground truth。

随后对 16 个 PDB pair 逐条完成机制分类，保留 {summary['curated_core_pairs']} 个 core benchmark pair；guide-dependent、lesion-specific、non-specific 和 transposase/substrate 样本被保留在 curated all table，但不进入 core specificity benchmark。基于 curated core 重新生成 benchmark_v0_2，共 {summary['benchmark_v0_2_rows']} 行，其中 single mutants {summary['single_mutants_v0_2']} 条、double mutants {summary['double_mutants_v0_2']} 条、random negatives {summary['random_negatives_v0_2']} 条。

分析方面，已重新输出 v0.2 figures，并完成 sequence-only proxy positional bias 检查，确认 k-mer/combined proxy 会产生位置效应，不能解释为生物学 specificity landscape。另建立了 Layer C experimental specificity pilot，从 JASPAR CORE 和 UniProt 获取 {summary['experimental_specificity_proteins']} 个 protein、{summary['experimental_specificity_rows']} 条 PFM-derived k-mer score 记录。下周重点是继续接入 raw PBM/HT-SELEX/CIS-BP 数据，并评估 protein-conditioned scoring interface。
"""
    (RESULTS_DIR / "WEEKLY_PROGRESS_V0_2.md").write_text(report, encoding="utf-8")


def main() -> None:
    curation = pd.read_csv(METADATA_DIR / "pdb_pair_curation.csv")
    pairs_v02 = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs_v0_2.csv")
    benchmark = pd.read_csv(PROCESSED_DIR / "benchmark_v0_2.csv", low_memory=False)
    experimental = pd.read_csv(PROCESSED_DIR / "experimental_specificity_small.csv")
    benchmark_counts = (
        benchmark.groupby("candidate_type")
        .size()
        .reset_index(name="count")
        .sort_values("candidate_type")
    )
    benchmark_counts.loc[len(benchmark_counts)] = {"candidate_type": "total", "count": int(len(benchmark))}
    summary = quality_summary(curation, pairs_v02, benchmark, experimental)
    pd.DataFrame([summary]).to_csv(RESULTS_TABLES / "benchmark_quality_summary_v0_2.csv", index=False)
    write_quality_report(summary, curation, benchmark_counts)
    write_weekly_progress(summary)
    print(f"Wrote {RESULTS_DIR / 'BENCHMARK_QUALITY_REPORT.md'}")
    print(f"Wrote {RESULTS_DIR / 'WEEKLY_PROGRESS_V0_2.md'}")
    print(pd.DataFrame([summary]).to_string(index=False))


if __name__ == "__main__":
    main()
