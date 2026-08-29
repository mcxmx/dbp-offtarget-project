from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.genome_candidates import download_genome, load_fasta_sequence, retrieve_candidates
from src.utils import ensure_dir, gc_content, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
RAW_GENOME_DIR = ensure_dir(ROOT / "data" / "raw" / "genome")
RESULTS_FIGURES = ensure_dir(ROOT / "results" / "figures")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")


def resolve_pair(pairs: pd.DataFrame, key: str) -> pd.Series:
    exact = pairs.loc[pairs["pair_id"] == key]
    if not exact.empty:
        return exact.iloc[0]
    by_pdb = pairs.loc[pairs["pdb_id"] == key]
    if not by_pdb.empty:
        return by_pdb.iloc[0]
    prefix = pairs.loc[pairs["pair_id"].astype(str).str.startswith(f"{key}_")]
    if not prefix.empty:
        return prefix.iloc[0]
    raise ValueError(f"Pair not found: {key}")


def write_weekly_progress() -> None:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")
    single = pd.read_csv(PROCESSED_DIR / "single_mutants.csv")
    double = pd.read_csv(PROCESSED_DIR / "double_mutants.csv")
    negatives = pd.read_csv(PROCESSED_DIR / "random_negatives.csv")
    benchmark = pd.read_csv(PROCESSED_DIR / "benchmark_v0_1.csv")
    genome = pd.read_csv(RESULTS_TABLES / "genome_candidates_demo.csv")
    lines = [
        "本周完成了一个可复现的 DBP-DNA off-target prototype。首先基于 RCSB PDB 整理了 "
        f"{len(pairs)} 条真实 protein-DNA 配对，并保留了结构条目、链 ID、PDB ID、来源 URL 和检索日期。",
        f"随后构建了 single mutant 数据集 {len(single)} 条、double mutant 数据集 {len(double)} 条，以及 GC-matched/random negatives 共 {len(negatives)} 条。",
        f"整合后的 benchmark_v0.1 共 {len(benchmark)} 条，已生成 sequence-only proxy baseline 和 preliminary figures。",
        f"genome scan demo 以 {CONFIG['genome_demo']['assembly']} {CONFIG['genome_demo']['chromosome']} 为目标，共返回 {len(genome)} 个候选位点。",
        "下周重点是扩充具有 quantitative specificity ground truth 的数据源，并把 protein-conditioned scoring 模型接口接到现有 benchmark 上。",
    ]
    text = "".join(lines)
    # Keep the report compact and factual.
    report = (
        "# Weekly Progress\n\n"
        "本周完成：\n"
        f"{lines[0]}\n"
        f"{lines[1]}\n"
        f"{lines[2]}\n"
        f"{lines[3]}\n\n"
        "下周计划：\n"
        f"{lines[4]}\n\n"
        "限制：当前结果仍是 sequence-only proxy，尚不能替代 protein-conditioned binding prediction。"
    )
    (ROOT / "results" / "WEEKLY_PROGRESS.md").write_text(report, encoding="utf-8-sig")


def main() -> None:
    genome_path = RAW_GENOME_DIR / "chr22.fa.gz"
    download_genome(CONFIG["genome_demo"]["fasta_url"], genome_path)
    genome = load_fasta_sequence(genome_path)
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")
    row = resolve_pair(pairs, CONFIG["genome_scan_pair_id"])
    target = normalize_sequence(row["target_dna"])
    candidates = retrieve_candidates(
        target,
        genome,
        int(CONFIG["genome_demo"]["seed_kmer_length"]),
        int(CONFIG["genome_demo"]["max_candidates"]),
        int(CONFIG["genome_demo"]["max_seed_hits_per_seed"]),
    )
    candidates.insert(0, "pair_id", row["pair_id"])
    candidates.insert(1, "target_dna", target)
    candidates.insert(2, "protein_name", row["protein_name"])
    candidates.insert(3, "protein_sequence", row["protein_sequence"])
    candidates["source"] = f"{CONFIG['genome_demo']['assembly']} {CONFIG['genome_demo']['chromosome']} scan"
    candidates["candidate_type"] = "genome_candidate"
    candidates.to_csv(RESULTS_TABLES / "genome_candidates_demo.csv", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "n_candidates", "value": int(len(candidates))},
            {"metric": "min_hamming", "value": int(candidates["hamming_distance"].min()) if not candidates.empty else None},
            {"metric": "median_hamming", "value": float(candidates["hamming_distance"].median()) if not candidates.empty else None},
            {"metric": "max_support_count", "value": int(candidates["support_count"].max()) if not candidates.empty else None},
        ]
    )
    summary.to_csv(RESULTS_TABLES / "genome_scan_summary.csv", index=False)

    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(8, 4))
    top = candidates.head(20).copy()
    if not top.empty:
        sns.barplot(data=top, x="rank", y="hamming_distance", hue="strand", dodge=False, ax=ax)
        ax.set_xlabel("rank")
        ax.set_ylabel("hamming distance")
        ax.set_title("Top genome scan candidates")
        fig.tight_layout()
        fig.savefig(RESULTS_FIGURES / "fig6_genome_scan_candidates.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    write_weekly_progress()
    print(f"Genome candidates written to {RESULTS_TABLES / 'genome_candidates_demo.csv'}")
    print(f"Weekly progress written to {ROOT / 'results' / 'WEEKLY_PROGRESS.md'}")


if __name__ == "__main__":
    main()
