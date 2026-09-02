from __future__ import annotations

import sys
from collections import defaultdict, deque
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.natural_pbm import EXPECTED_8MER_RC_CLASSES
from src.utils import ensure_dir, project_root


ROOT = project_root()
INTERIM = ROOT / "data" / "interim" / "v0_4_1"
PROCESSED = ensure_dir(ROOT / "data" / "processed" / "v0_4_1")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4_1")
RESULTS = ensure_dir(ROOT / "results" / "v0_4_1")
TABLES = ensure_dir(RESULTS / "tables")
FIGURES = ensure_dir(RESULTS / "figures")
DOCS = ensure_dir(ROOT / "docs" / "v0_4_1")
TODAY = date.today().isoformat()
RNG = np.random.default_rng(42)


FUSION_PATTERNS = ("GST-", "MAML1", "NOTCH")
HOMEODOMAIN = {"Abd-B", "Bap", "Eve", "Lbl", "Msh", "Ptx1", "Six4", "Slou", "Tin", "Ubx", "Lmd", "CEH-22", "Ceh-22"}
BZIP = {"GCN4", "Jun_Fos", "Cad1", "Cin5", "Hac1", "Yap3", "Cst6", "Sko1"}
ZINC_CLUSTER = {"Ecm22", "Hap1", "Pdr3", "Upc2", "Sut1"}
C2H2 = {"BCL11A", "BCL11B", "Zif268", "Zap1", "Mot3", "Stb4", "Stb5", "Vhr1", "Vhr2"}
BHLH = {"HLH-1", "Cbf1"}
MYB = {"Rap1", "Abf1"}


def coarse_family(name: str) -> str:
    if name in HOMEODOMAIN:
        return "homeodomain"
    if name in BZIP:
        return "bZIP_or_AP1_like"
    if name in ZINC_CLUSTER:
        return "fungal_Zn2Cys6_like"
    if name in C2H2:
        return "C2H2_zinc_finger_like"
    if name in BHLH:
        return "bHLH"
    if name in MYB:
        return "Myb_like"
    return "unknown_or_publication_specific"


def is_fusion_or_unclear(name: str) -> bool:
    return any(pattern in name for pattern in FUSION_PATTERNS) or "contig8mers" in name


def build_quality(long_df: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.set_index("natural_protein_id", drop=False)
    rows = []
    for protein_id, group in long_df.groupby("natural_protein_id", sort=True):
        row = meta.loc[protein_id]
        name = str(group["protein_name"].iloc[0])
        n_rc = int(group["canonical_rc"].nunique())
        complete = n_rc == EXPECTED_8MER_RC_CLASSES
        valid = bool(group["dna_valid"].all()) and bool(group["rc_source_matches"].all())
        has_seq = pd.notna(row.get("protein_sequence"))
        unclear = is_fusion_or_unclear(name)
        conf = str(row.get("sequence_confidence", "low"))
        if not valid or unclear:
            quality = "exclude"
            reason = "invalid DNA/RC check or complex/fusion/unclear factor label"
        elif not complete:
            quality = "low_confidence"
            reason = "incomplete contiguous 8-mer RC-class coverage"
        elif not has_seq:
            quality = "low_confidence"
            reason = "complete PBM profile but no conservative protein sequence recovery"
        elif conf == "high":
            quality = "high_confidence"
            reason = "complete PBM profile and high-confidence reference protein sequence"
        else:
            quality = "usable"
            reason = "complete PBM profile and medium-confidence reference protein sequence"
        rows.append(
            {
                "natural_protein_id": protein_id,
                "protein_name": name,
                "n_rows": int(len(group)),
                "n_unique_8mers": int(group["dna_sequence"].nunique()),
                "n_rc_classes": n_rc,
                "n_missing_rc_classes": int(EXPECTED_8MER_RC_CLASSES - n_rc),
                "n_duplicate_rows": int(len(group) - group["dna_sequence"].nunique()),
                "n_experiments": int(group["experiment_id"].nunique()),
                "n_replicate_groups": int(group["replicate_id"].nunique()),
                "dna_qc_pass": bool(group["dna_valid"].all()),
                "rc_qc_pass": bool(group["rc_source_matches"].all()),
                "has_protein_sequence": bool(has_seq),
                "sequence_confidence": conf,
                "quality_level": quality,
                "quality_reason": reason,
                "training_candidate": quality in {"high_confidence", "usable"},
            }
        )
    qc = pd.DataFrame(rows)
    qc.to_csv(TABLES / "natural_pbm_qc_summary.csv", index=False)
    return qc


def build_benchmark(long_df: pd.DataFrame, meta: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()
    meta["protein_family"] = meta["protein_name"].map(coarse_family)
    manifest = pd.read_csv(METADATA / "natural_pbm_files.csv")
    url_map = dict(zip(manifest["dataset_code"], manifest["url"]))
    meta["dataset_code"] = meta["natural_protein_id"].astype(str).str.split(":").str[1]
    meta["source_url"] = meta["dataset_code"].map(url_map).fillna(meta.get("source_url", pd.NA))
    meta.to_csv(METADATA / "natural_pbm_proteins.csv", index=False)
    base_path = PROCESSED / "natural_pbm_benchmark_v0_4_1.parquet"
    if not base_path.exists():
        raise RuntimeError("Initial consensus benchmark is missing; run 01_build_natural_pbm_benchmark.py first.")
    base = pd.read_parquet(base_path)
    keep_cols = [
        "natural_protein_id",
        "protein_name",
        "protein_sequence",
        "sequence_length",
        "sequence_type",
        "sequence_match_to_assay",
        "uniprot_id",
        "species",
        "protein_family",
    ]
    meta_small = meta.loc[:, [col for col in keep_cols if col in meta.columns]].rename(
        columns={"natural_protein_id": "protein_id"}
    )
    base = base.drop(
        columns=[
            col
            for col in [
                "protein_sequence",
                "sequence_length",
                "sequence_type",
                "sequence_match_to_assay",
                "uniprot_id",
                "species",
                "protein_family",
                "quality_level",
                "training_candidate",
            ]
            if col in base.columns
        ],
        errors="ignore",
    )
    merged = base.merge(meta_small, on=["protein_id", "protein_name"], how="left")
    merged = merged.merge(
        qc[["natural_protein_id", "quality_level", "training_candidate"]].rename(
            columns={"natural_protein_id": "protein_id"}
        ),
        on="protein_id",
        how="left",
    )
    merged = merged[merged["training_candidate"]].copy()
    merged["experimental_percentile"] = merged.groupby("protein_id")["experimental_score"].rank(pct=True, ascending=True)
    merged["experimental_rank"] = merged.groupby("protein_id")["experimental_score"].rank(method="average", ascending=False)
    merged["within_protein_normalized_rank"] = 1.0 - (
        (merged["experimental_rank"] - 1.0)
        / merged.groupby("protein_id")["experimental_rank"].transform("max").clip(lower=1.0)
    )
    out = pd.DataFrame(
        {
            "protein_id": merged["protein_id"],
            "protein_name": merged["protein_name"],
            "protein_sequence": merged["protein_sequence"],
            "sequence_length": merged["sequence_length"],
            "sequence_type": merged["sequence_type"],
            "sequence_match_to_assay": merged["sequence_match_to_assay"],
            "uniprot_id": merged["uniprot_id"],
            "protein_family": merged["protein_family"],
            "species": merged["species"],
            "dna_sequence": merged["canonical_rc"],
            "dna_length": merged["dna_length"],
            "canonical_rc": merged["canonical_rc"],
            "experimental_score": merged["experimental_score"],
            "experimental_percentile": merged["experimental_percentile"],
            "experimental_rank": merged["experimental_rank"],
            "within_protein_normalized_rank": merged["within_protein_normalized_rank"],
            "experiment_id": merged["experiment_id"],
            "assay_type": merged["assay_type"],
            "score_type": merged["score_type"] if "score_type" in merged.columns else merged["experimental_score_type"],
            "quality_level": merged["quality_level"],
            "source": merged["source"],
            "source_file": merged["source_file"] if "source_file" in merged.columns else pd.NA,
            "retrieval_date": TODAY,
        }
    )
    out.to_parquet(PROCESSED / "natural_pbm_benchmark_v0_4_1.parquet", index=False)
    out.head(5000).to_csv(TABLES / "natural_pbm_benchmark_v0_4_1_preview.csv", index=False)
    return out


def sequence_identity(seq_a: str, seq_b: str) -> float:
    if seq_a == seq_b:
        return 1.0
    k = 3
    if len(seq_a) < k or len(seq_b) < k:
        return 0.0
    a = {seq_a[i : i + k] for i in range(len(seq_a) - k + 1)}
    b = {seq_b[i : i + k] for i in range(len(seq_b) - k + 1)}
    if not a or not b:
        return 0.0
    # This fast proxy is used only for v0.4.1 split hygiene. A later larger
    # release should replace it with MMseqs2/CD-HIT on Linux.
    return float(len(a & b) / len(a | b))


def connected_components(ids: list[str], matrix: pd.DataFrame, threshold: float) -> dict[str, str]:
    graph = defaultdict(set)
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            if float(matrix.loc[a, b]) >= threshold:
                graph[a].add(b)
                graph[b].add(a)
    seen = set()
    out = {}
    cluster_idx = 0
    for start in ids:
        if start in seen:
            continue
        cluster_idx += 1
        queue = deque([start])
        members = []
        seen.add(start)
        while queue:
            item = queue.popleft()
            members.append(item)
            for nxt in graph[item]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        label = f"cluster_{threshold:.2f}_{cluster_idx:03d}"
        for member in members:
            out[member] = label
    return out


def cluster_and_split(bench: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    proteins = (
        bench[
            [
                "protein_id",
                "protein_name",
                "protein_sequence",
                "protein_family",
                "species",
                "quality_level",
            ]
        ]
        .drop_duplicates("protein_id")
        .sort_values("protein_id")
        .reset_index(drop=True)
    )
    ids = proteins["protein_id"].tolist()
    matrix = pd.DataFrame(np.eye(len(ids)), index=ids, columns=ids)
    seqs = dict(zip(proteins["protein_id"], proteins["protein_sequence"]))
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            ident = sequence_identity(str(seqs[a]), str(seqs[b]))
            matrix.loc[a, b] = ident
            matrix.loc[b, a] = ident
    matrix.to_csv(TABLES / "natural_protein_sequence_identity_matrix.csv")
    clusters = proteins.copy()
    for pct in [0.30, 0.40, 0.50, 0.60]:
        clusters[f"cluster_{int(pct * 100)}"] = clusters["protein_id"].map(connected_components(ids, matrix, pct))
    clusters.to_csv(METADATA / "natural_protein_clusters.csv", index=False)

    cluster_labels = sorted(clusters["cluster_40"].unique())
    RNG.shuffle(cluster_labels)
    cluster_sizes = clusters.groupby("cluster_40")["protein_id"].nunique().to_dict()
    totals = {"train": 0, "validation": 0, "natural_test": 0}
    targets = {"train": 0.70 * len(ids), "validation": 0.15 * len(ids), "natural_test": 0.15 * len(ids)}
    assignment = {}
    for cluster in sorted(cluster_labels, key=lambda c: cluster_sizes[c], reverse=True):
        split = min(totals, key=lambda s: totals[s] / max(targets[s], 1))
        assignment[cluster] = split
        totals[split] += cluster_sizes[cluster]
    splits = clusters[["protein_id", "protein_name", "protein_family", "species", "cluster_40"]].copy()
    splits["split"] = splits["cluster_40"].map(assignment)
    counts = (
        splits.groupby("split")
        .agg(n_proteins=("protein_id", "nunique"), n_families=("protein_family", "nunique"))
        .reset_index()
    )
    dna_counts = bench.merge(splits[["protein_id", "split"]], on="protein_id", how="left").groupby("split").size()
    counts["n_protein_dna_units"] = counts["split"].map(dna_counts).astype(int)
    splits.to_csv(METADATA / "natural_pbm_splits.csv", index=False)
    counts.to_csv(TABLES / "natural_pbm_split_summary.csv", index=False)
    return clusters, splits, counts


def write_docs(qc: pd.DataFrame, bench: pd.DataFrame, splits: pd.DataFrame, split_counts: pd.DataFrame) -> None:
    DOCS.joinpath("NATURAL_PBM_QC_RULES.md").write_text(
        f"""# v0.4.1 Natural PBM QC Rules

Audit date: {TODAY}

## Unit Definitions

- Natural PBM source rows are UniPROBE processed contiguous 8-mer E-score rows.
- The independent DNA unit used for splitting/evaluation is the reverse-complement canonical 8-mer class.
- A complete contiguous 8-mer profile is expected to contain {EXPECTED_8MER_RC_CLASSES:,} RC classes.

## Quality Levels

- `high_confidence`: complete RC-class coverage, DNA/RC QC pass, clear non-fusion protein label, and high-confidence UniProt reference sequence.
- `usable`: complete RC-class coverage, DNA/RC QC pass, clear non-fusion protein label, and medium-confidence UniProt reference sequence.
- `low_confidence`: complete/incomplete PBM profile that lacks a conservative protein sequence, or incomplete 8-mer coverage.
- `exclude`: invalid DNA/RC QC or complex/fusion/unclear labels.

Only `high_confidence` and `usable` rows enter `natural_pbm_benchmark_v0_4_1.parquet`.
Full-length UniProt sequences are retained with `sequence_match_to_assay=false` unless construct-level sequence evidence is available.

Protein clusters in v0.4.1 are generated from a fast amino-acid 3-mer Jaccard proxy to keep obvious duplicate and near-duplicate reference sequences in the same split. This is a split-hygiene proxy, not a replacement for a future MMseqs2/CD-HIT homology audit.
""",
        encoding="utf-8",
    )
    report = f"""# v0.4.1 Natural PBM Benchmark Report

Audit date: {TODAY}

## Summary

- Source: UniPROBE processed contiguous 8-mer E-score profiles.
- Final train/evaluation benchmark proteins: {bench['protein_id'].nunique():,}
- Protein families/coarse classes: {bench[['protein_id', 'protein_family']].drop_duplicates()['protein_family'].nunique():,}
- Species represented: {bench[['protein_id', 'species']].drop_duplicates()['species'].nunique():,}
- Protein-DNA units: {len(bench):,} protein-RC-class rows.
- DNA unit: 8-mer reverse-complement equivalence class.
- Score: UniPROBE contiguous 8-mer E-score; used for per-protein ranking, not cross-protein absolute affinity.
- Protein sequence completeness in final benchmark: {bench[['protein_id', 'protein_sequence']].drop_duplicates()['protein_sequence'].notna().mean():.1%}
- Construct sequence completeness: 0.0%; reference full-length sequences are used and flagged as not assay-construct matched.

## QC Counts

{qc['quality_level'].value_counts().to_string()}

## Cluster-Aware Split

{split_counts.to_string(index=False)}

## Training Readiness

The benchmark is sufficient for a first simple protein-conditioned baseline because it has >50 proteins with complete 8-mer profiles and conservative reference sequences. It is not yet sufficient for claims about assay-free biological OOD because natural UniPROBE 8-mer profiles and designed GSE237017 uPBM 7-mer profiles differ in k-mer length, protocol, and score processing.
"""
    (RESULTS / "NATURAL_PBM_BENCHMARK_REPORT.md").write_text(report, encoding="utf-8")


def make_figures(qc: pd.DataFrame, bench: pd.DataFrame, clusters: pd.DataFrame, rep: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")

    per_source = bench[["protein_id", "source"]].drop_duplicates()["source"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    per_source.plot(kind="barh", ax=ax, color="#3A6EA5")
    ax.set_xlabel("Proteins in final benchmark")
    ax.set_ylabel("UniPROBE publication")
    ax.set_title("Natural PBM Dataset Overview")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_v0_4_1_1_natural_pbm_dataset_overview.png", dpi=300)
    plt.close(fig)

    family_counts = clusters["protein_family"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(7, 4))
    family_counts.plot(kind="barh", ax=ax, color="#6B8E23")
    ax.set_xlabel("Proteins")
    ax.set_ylabel("Coarse family annotation")
    ax.set_title("Natural PBM Protein Family Distribution")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_v0_4_1_2_natural_pbm_family_distribution.png", dpi=300)
    plt.close(fig)

    rep_sub = rep.dropna(subset=["spearman"]).copy()
    fig, ax = plt.subplots(figsize=(5.5, 4))
    if not rep_sub.empty:
        ax.scatter(rep_sub["pearson"], rep_sub["spearman"], s=45, color="#7A4E9D")
    ax.set_xlabel("Pearson replicate correlation")
    ax.set_ylabel("Spearman replicate correlation")
    ax.set_title("Natural PBM Replicate QC")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_v0_4_1_3_natural_replicate_qc.png", dpi=300)
    plt.close(fig)


def main() -> None:
    long_df = pd.read_parquet(INTERIM / "natural_pbm_long.parquet")
    meta = pd.read_csv(METADATA / "natural_pbm_proteins.csv")
    rep = pd.read_csv(TABLES / "natural_pbm_replicate_qc.csv")
    qc = build_quality(long_df, meta)
    bench = build_benchmark(long_df, meta, qc)
    clusters, splits, split_counts = cluster_and_split(bench)
    write_docs(qc, bench, splits, split_counts)
    make_figures(qc, bench, clusters, rep)
    print(f"final natural proteins={bench['protein_id'].nunique()} units={len(bench)}")
    print(qc["quality_level"].value_counts().to_string())
    print(split_counts.to_string(index=False))


if __name__ == "__main__":
    main()
