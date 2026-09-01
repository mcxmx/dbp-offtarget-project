from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
DOCS_DIR = ensure_dir(ROOT / "docs" / "v0_3")
RESULTS_DIR = ensure_dir(ROOT / "results" / "v0_3")
TABLES_DIR = ensure_dir(RESULTS_DIR / "tables")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed" / "v0_3")
METADATA_DIR = ensure_dir(ROOT / "metadata" / "v0_3")


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


def short_float(value: float | str, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def load_tables() -> dict[str, pd.DataFrame]:
    return {
        "samples": pd.read_csv(METADATA_DIR / "gse237017_samples.csv", dtype={"replicate": str}),
        "manifest": pd.read_csv(METADATA_DIR / "gse237017_file_manifest.csv"),
        "sequences": pd.read_csv(METADATA_DIR / "designed_dbp_sequences.csv"),
        "targets": pd.read_csv(METADATA_DIR / "designed_dbp_targets.csv"),
        "coverage": pd.read_csv(TABLES_DIR / "sample_coverage_qc.csv", dtype={"replicate": str}),
        "replicate": pd.read_csv(TABLES_DIR / "replicate_qc.csv"),
        "rc": pd.read_csv(TABLES_DIR / "reverse_complement_qc.csv"),
        "specificity": pd.read_csv(TABLES_DIR / "protein_specificity_summary.csv"),
        "target_rank": pd.read_csv(TABLES_DIR / "target_rank_summary.csv"),
        "baseline": pd.read_csv(TABLES_DIR / "designed_dbp_sequence_baseline.csv"),
        "failures": pd.read_csv(TABLES_DIR / "sequence_similarity_failure_candidates.csv"),
        "benchmark": pd.read_parquet(PROCESSED_DIR / "designed_dbp_upbm_v0_3.parquet"),
    }


def summarize(tables: dict[str, pd.DataFrame]) -> dict[str, object]:
    samples = tables["samples"]
    replicate = tables["replicate"]
    benchmark = tables["benchmark"]
    valid_reps = replicate[replicate["qc_status"] == "replicate_pair"].copy()
    e_rep = valid_reps[valid_reps["score_type"] == "e_score"].copy()
    single_rep = sorted(replicate.loc[replicate["qc_status"] == "single_replicate_only", "protein_id"].unique())
    return {
        "n_samples": int(samples["gsm_id"].nunique()),
        "n_proteins": int(samples["protein_id"].nunique()),
        "n_replicate_pairs": int(len(e_rep)),
        "replicate_proteins": "; ".join(sorted(e_rep["protein_id"].unique())),
        "single_replicate_proteins": "; ".join(single_rep),
        "n_measurements": int(len(benchmark)),
        "n_unique_7mers_per_protein": int(benchmark.groupby("protein_id")["dna_7mer"].nunique().min()),
        "protein_sequences_recovered": int(tables["sequences"]["protein_sequence"].notna().sum()),
        "targets_recovered": int(tables["targets"]["intended_target_dna"].notna().sum()),
        "e_score_pearson_min": float(e_rep["pearson_correlation"].min()),
        "e_score_pearson_median": float(e_rep["pearson_correlation"].median()),
        "e_score_pearson_max": float(e_rep["pearson_correlation"].max()),
        "e_score_spearman_min": float(e_rep["spearman_correlation"].min()),
        "e_score_spearman_median": float(e_rep["spearman_correlation"].median()),
        "e_score_spearman_max": float(e_rep["spearman_correlation"].max()),
        "baseline_spearman_min": float(tables["baseline"]["spearman"].min()),
        "baseline_spearman_median": float(tables["baseline"]["spearman"].median()),
        "baseline_spearman_max": float(tables["baseline"]["spearman"].max()),
        "failure_candidates": int(len(tables["failures"])),
    }


def write_score_definitions() -> None:
    text = """# PBM Score Definitions for GSE237017

Source: GEO accession GSE237017, sample-level processed uPBM files, and GEO data processing text.

## Experiment

GSE237017 contains universal protein-binding microarray experiments for designed DNA-binding proteins DBP1, DBP3, DBP5, DBP6, DBP9, DBP35, and DBP48. GEO describes the arrays as 15K Agilent dsDNA arrays designed so that all possible 9-bp sequences are covered, with every 7-mer represented in at least 16 spots.

## Data Processing Recorded by GEO

The GEO sample metadata states that Alexa 488 signals were position-adjusted to correct microarray non-uniformity. A Seed-and-wobble algorithm was then used to compute enrichment scores and median intensities for all possible 7-mers. GEO defines the processed values as E-scores, median intensities, and z-scores for 7-mers.

## E-score

Recommended v0.3 primary score: `PBM E-score`, stored as `experimental_score_primary`.

Interpretation:

- Experimental PBM specificity/enrichment score for 7-mer ranking within a protein.
- Useful for per-protein specificity landscape analysis.
- Suitable as the first primary target for designed-DBP external benchmark ranking.

Not justified:

- Not Kd.
- Not binding free energy.
- Not binding probability.
- Not directly comparable as absolute affinity across different proteins.

## Median Intensity

Interpretation:

- Processed microarray signal intensity summary for 7-mers.
- Useful as a secondary assay-derived measurement.

Not justified:

- Not directly calibrated affinity.
- Not automatically comparable across proteins or experiments without normalization.

## Z-score

Interpretation:

- Standardized processed score included by GEO for 7-mers.
- Useful as a secondary consistency/QC measure.

Not justified:

- Not a probability or Kd.
- Not used as the primary v0.3 score unless future analysis shows it is preferable.

## Reverse Complement Handling

The processed 7-mer files contain two 7-mer columns. In v0.3 the parser explicitly expands both the primary and reverse-complement companion columns, assigning the same GEO-provided scores to both 7-mer orientations. `reverse_complement_qc.csv` confirms zero difference between paired reverse-complement scores after expansion.
"""
    (DOCS_DIR / "PBM_SCORE_DEFINITIONS.md").write_text(text, encoding="utf-8")


def write_dataset_card(tables: dict[str, pd.DataFrame], stats: dict[str, object]) -> None:
    sample_table = tables["samples"][["gsm_id", "protein_id", "protein_concentration", "replicate", "sample_title"]]
    target_table = tables["targets"][["protein_id", "intended_target_dna", "target_length", "target_id", "target_context", "confidence"]]
    text = f"""# Designed DBP uPBM Dataset Card v0.3

## Dataset Name

Designed DBP uPBM experimental specificity benchmark v0.3.

## Source

- GEO accession: GSE237017
- GEO title: Computational design of sequence-specific DNA-binding proteins
- Paper reference: DOI 10.1038/s41594-025-01669-4; PMID 40940539
- Primary GEO URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237017

## Proteins

The benchmark contains {stats['n_proteins']} designed DNA-binding proteins: DBP1, DBP3, DBP5, DBP6, DBP9, DBP35, and DBP48.

Protein sequences recovered from official supplementary material: {stats['protein_sequences_recovered']}/7.
Intended target DNA sequences recovered from official supplementary material: {stats['targets_recovered']}/7.

## Samples and Replicates

Usable GEO samples: {stats['n_samples']}.
Replicate pairs with E-score QC: {stats['n_replicate_pairs']}.
Proteins with replicate samples: {stats['replicate_proteins']}.
Single-replicate-only proteins: {stats['single_replicate_proteins']}.

{markdown_table(sample_table)}

## Experimental Protocol

Universal protein-binding microarray experiments were performed for small recombinant designed DNA-binding proteins. The processed files provide 7-mer E-scores, median intensities, and z-scores.

## DNA Sequence Space

Each sample covers {stats['n_unique_7mers_per_protein']} unique 7-mers after explicit expansion of the reverse-complement companion column. The final consensus benchmark contains {stats['n_measurements']} protein-7mer measurements.

## Score Definitions

Primary score: PBM E-score, stored as `experimental_score_primary`.

Secondary scores: `median_intensity_mean`, `z_score_mean`.

The primary score is an experimental PBM specificity/enrichment score for per-protein ranking. It is not a Kd, binding probability, binding free energy, or absolute cross-protein affinity.

## Processing Steps

1. Downloaded GSE237017 family SOFT metadata and GSM supplementary files.
2. Parsed sample title, DBP ID, concentration, replicate, platform, and supplementary URLs from GEO metadata.
3. Downloaded all 12 processed 7-mer files and 12 raw spot-data files.
4. Preserved SHA256 and file size for every raw file.
5. Parsed processed 7-mer tables and expanded primary/reverse-complement 7-mer columns.
6. Built replicate-level consensus by mean E-score within each protein and 7-mer.
7. Added per-protein rank, percentile, and within-protein z-score without overwriting raw E-score.

## QC

- E-score replicate Pearson range: {short_float(stats['e_score_pearson_min'])} to {short_float(stats['e_score_pearson_max'])}; median {short_float(stats['e_score_pearson_median'])}.
- E-score replicate Spearman range: {short_float(stats['e_score_spearman_min'])} to {short_float(stats['e_score_spearman_max'])}; median {short_float(stats['e_score_spearman_median'])}.
- Reverse-complement score differences are zero after expansion.
- Sample coverage is complete for all processed samples after RC-column expansion.

## Intended Targets

{markdown_table(target_table)}

## Known Limitations

- Processed specificity is 7-mer based; full target DNA binding is not directly measured by these tables.
- DBP5 and DBP48 have one GEO sample each in the parsed metadata, so replicate QC cannot be computed for them.
- Absolute E-scores should not be compared across proteins as binding affinity.
- The benchmark is in vitro uPBM, not in vivo genomic binding.

## Recommended Use

- Designed DBP specificity ranking.
- External validation for protein-conditioned DNA specificity models.
- Out-of-distribution evaluation of models trained on natural DBPs.

## Prohibited Interpretation

- Absolute affinity regression across proteins.
- Direct Kd or free-energy interpretation.
- Direct in vivo genomic binding claim.
- Claiming full-target affinity from overlapping 7-mer summaries.
"""
    (DOCS_DIR / "DESIGNED_DBP_UPBM_DATASET_CARD.md").write_text(text, encoding="utf-8")


def write_model_split_plan() -> None:
    text = """# Model Split Plan for v0.4+

GSE237017 designed DBPs should be protected as an external designed-protein benchmark rather than mixed into natural transcription-factor training.

## Split A: Natural Protein Random Split

Use only natural DBP experimental specificity datasets. Randomly split proteins after deduplicating identical protein sequences and highly similar motifs. This split is useful for debugging but should not be the headline generalization result.

## Split B: Protein-Family Split

Cluster natural proteins by sequence identity, domain annotation, and family labels where available. Hold out entire protein families to measure family-level generalization.

## Split C: Designed DBP Zero-Shot External Test

Train on natural DBPs only. Evaluate DBP1, DBP3, DBP5, DBP6, DBP9, DBP35, and DBP48 from GSE237017 as external/OOD designed DBPs. This is one of the key paper-level experiments: testing whether natural-trained protein-DNA models suffer an OOD performance drop on de novo designed binders.

## Split D: Leave-One-Designed-DBP-Out

After a zero-shot evaluation is established, optionally train/calibrate on six designed DBPs and test on the held-out designed DBP. This evaluates transfer within the designed-binder distribution, but it should not replace Split C.

## Leakage Controls

- Do not allow identical protein sequences across train/test.
- Cluster related natural proteins by family.
- Track DNA target similarity and motif similarity separately.
- Keep GSE237017 raw 7-mer scores separate from synthetic mutation benchmarks.
- Report per-protein metrics rather than pooling absolute scores across proteins.

## Recommended Metrics

- Per-protein Spearman correlation.
- Top-k enrichment within protein.
- Ranking quality for high-scoring 7-mers.
- OOD drop from natural-family held-out sets to designed DBPs.
"""
    (DOCS_DIR / "MODEL_SPLIT_PLAN.md").write_text(text, encoding="utf-8")


def write_progress_and_go_no_go(tables: dict[str, pd.DataFrame], stats: dict[str, object]) -> None:
    target_rank = tables["target_rank"].copy()
    baseline = tables["baseline"].copy()
    target_rank["best_target_percentile_label"] = target_rank["best_target_percentile"].map(lambda x: short_float(x, 3))
    target_table = target_rank[
        ["protein_id", "best_target_7mer", "best_target_percentile_label", "mean_target_7mer_escore"]
    ].rename(columns={"best_target_percentile_label": "best_target_percentile"})
    baseline_summary = baseline.groupby("protein_id")["spearman"].agg(["min", "median", "max"]).reset_index()
    dbp48_percentile = target_rank.loc[target_rank["protein_id"] == "DBP48", "best_target_percentile_label"].iloc[0]

    progress = f"""# v0.3 Progress Report

本轮建立了 designed DBP experimental specificity benchmark v0.3，核心数据源为 GEO GSE237017。已程序化解析 series/GSM metadata，并下载 12 个 usable uPBM samples，覆盖 7 个 designed DBP：DBP1、DBP3、DBP5、DBP6、DBP9、DBP35、DBP48。每个 sample 的 processed 7-mer 文件原始为 8192 行、两个 7-mer 列；解析时显式展开 reverse-complement companion column 后，每个 sample 覆盖 16,384 个 unique 7-mers，缺失数为 0。最终 consensus benchmark 包含 {stats['n_measurements']} 条 protein-7mer experimental measurements。

Replicate QC 显示 DBP1、DBP3、DBP6、DBP9、DBP35 有 replicate；DBP5 和 DBP48 为 single replicate only。E-score replicate Pearson 范围为 {short_float(stats['e_score_pearson_min'])}-{short_float(stats['e_score_pearson_max'])}，中位数 {short_float(stats['e_score_pearson_median'])}；Spearman 范围为 {short_float(stats['e_score_spearman_min'])}-{short_float(stats['e_score_spearman_max'])}，中位数 {short_float(stats['e_score_spearman_median'])}。7 个 DBP 的 protein sequence 和 intended target sequence 均已从官方 supplementary workbook 恢复，confidence 为 high。

Target rank 分析只基于 intended target 的 overlapping 7-mers，不解释为 full-target affinity。多数设计的 best target-derived 7-mer 位于较高 percentile，但 DBP48 的 best target-derived 7-mer percentile 为 {dbp48_percentile}，相对较弱。Sequence-only baseline 与 PBM E-score 的 per-protein Spearman 整体较低，中位数接近 0；同时发现 {stats['failure_candidates']} 个 top 1% E-score 但 target-derived 7-mer similarity 不高的候选，说明仅靠 DNA sequence similarity 难以解释 designed DBP specificity landscape。

当前最大限制是 PBM processed score 是 7-mer 级别，不能直接代表完整 target DNA 的 binding affinity；DBP5/DBP48 缺少 replicate；E-score 也不能跨 protein 当作绝对 affinity。下一步应进入 v0.4 protein-conditioned baseline：以 protein sequence、intended target-derived context 和 candidate 7-mer 为输入，先做非神经或轻量模型的 per-protein ranking baseline，再和 sequence-only baseline 比较。
"""
    (RESULTS_DIR / "V0_3_PROGRESS.md").write_text(progress, encoding="utf-8")

    go = f"""# GO / NO-GO for v0.3 Designed DBP uPBM Benchmark

Decision: CONDITIONAL GO

## 1. Is GSE237017 sufficient as a designed-DBP external benchmark?

Yes, conditionally. It provides {stats['n_samples']} usable uPBM samples across {stats['n_proteins']} designed DBPs and {stats['n_measurements']} protein-7mer measurements with GEO provenance. It is suitable for external/OOD ranking evaluation, not absolute affinity regression.

## 2. Are replicates sufficiently consistent?

Partially. E-score replicate Pearson median is {short_float(stats['e_score_pearson_median'])} and Spearman median is {short_float(stats['e_score_spearman_median'])}. This is usable for ranking analyses with QC caveats. DBP5 and DBP48 have single replicate only in GEO metadata.

## 3. Were all designed DBP protein sequences recovered?

Yes. {stats['protein_sequences_recovered']}/7 protein sequences were recovered from the official Nature supplementary workbook.

## 4. Were intended targets recovered?

Yes. {stats['targets_recovered']}/7 intended target DNA sequences were recovered from the official Nature supplementary workbook.

## 5. How does sequence-only similarity correlate with experimental specificity?

Weakly. Per-protein sequence-only Spearman correlations range from {short_float(stats['baseline_spearman_min'])} to {short_float(stats['baseline_spearman_max'])}, with median {short_float(stats['baseline_spearman_median'])}. This supports using v0.3 to test protein-conditioned models.

## 6. Are there high-score 7-mers not explained by simple sequence similarity?

Yes. The v0.3 analysis found {stats['failure_candidates']} candidate rows with top 1% PBM E-score and hamming similarity to target-derived 7-mers at or below the protein median.

## 7. Is the project ready for protein-conditioned baselines?

Conditionally yes. The dataset is ready for per-protein ranking baselines and zero-shot designed-DBP evaluation. It is not ready for calibrated affinity or uncertainty claims.

## Main Limitation

The current benchmark is 7-mer in vitro uPBM specificity data. Full intended target sequences are longer than 7 bp, so target-rank summaries use overlapping 7-mers and must not be treated as full-target binding affinity.

## Recommended Next Baseline

Run a protein-conditioned but non-neural baseline first: encode the designed DBP sequence and candidate 7-mer, evaluate per-protein Spearman and top-k enrichment, and compare against the v0.3 sequence-only baseline.
"""
    (RESULTS_DIR / "GO_NO_GO.md").write_text(go, encoding="utf-8")

    pd.DataFrame([stats]).to_csv(TABLES_DIR / "v0_3_dataset_summary.csv", index=False)
    target_table.to_csv(TABLES_DIR / "v0_3_target_rank_compact.csv", index=False)
    baseline_summary.to_csv(TABLES_DIR / "v0_3_sequence_baseline_compact.csv", index=False)


def main() -> None:
    tables = load_tables()
    stats = summarize(tables)
    write_score_definitions()
    write_dataset_card(tables, stats)
    write_model_split_plan()
    write_progress_and_go_no_go(tables, stats)
    print(pd.DataFrame([stats]).to_string(index=False))


if __name__ == "__main__":
    main()
