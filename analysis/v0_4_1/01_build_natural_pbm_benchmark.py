from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.natural_pbm import (
    EXPECTED_8MER_RC_CLASSES,
    consensus_by_experiment,
    parse_uniprobe_archives,
    reverse_complement,
)
from src.utils import ensure_dir, project_root


ROOT = project_root()
RAW = ROOT / "data" / "raw" / "v0_4_1"
INTERIM = ensure_dir(ROOT / "data" / "interim" / "v0_4_1")
PROCESSED = ensure_dir(ROOT / "data" / "processed" / "v0_4_1")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4_1")
RESULTS = ensure_dir(ROOT / "results" / "v0_4_1")
TABLES = ensure_dir(RESULTS / "tables")
FIGURES = ensure_dir(RESULTS / "figures")
TODAY = date.today().isoformat()


def build_metadata_frame() -> pd.DataFrame:
    manifest = pd.read_csv(METADATA / "natural_pbm_files.csv")
    long_df = parse_uniprobe_archives(manifest)
    if long_df.empty:
        raise RuntimeError("No UniPROBE natural PBM rows could be parsed")
    long_path = INTERIM / "natural_pbm_long.parquet"
    long_df.to_parquet(long_path, index=False)
    return long_df


def protein_metadata(long_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein_id, group in long_df.groupby("natural_protein_id", sort=True):
        rows.append(
            {
                "natural_protein_id": protein_id,
                "protein_name": group["protein_name"].iloc[0],
                "protein_sequence": pd.NA,
                "sequence_length": pd.NA,
                "sequence_type": "unknown",
                "sequence_match_to_assay": False,
                "uniprot_id": pd.NA,
                "species": pd.NA,
                "protein_family": pd.NA,
                "experiment_id": ";".join(sorted(group["experiment_id"].unique())),
                "publication": ";".join(sorted(group["publication"].unique())),
                "doi": pd.NA,
                "database": "UniPROBE",
                "assay_type": "universal protein-binding microarray",
                "kmer_length": 8,
                "score_type": "UniPROBE contiguous 8-mer E-score",
                "replicate_group": "multi" if group["replicate_id"].nunique() > 1 else "single_measurement",
                "source_url": pd.NA,
                "retrieval_date": TODAY,
                "notes": "Protein sequence is not yet recovered from this publication-level 8-mer archive; construct-level provenance to be added separately.",
            }
        )
    meta = pd.DataFrame(rows)
    meta.to_csv(METADATA / "natural_pbm_proteins.csv", index=False)
    return meta


def quality_tables(long_df: pd.DataFrame) -> pd.DataFrame:
    qc_rows = []
    for protein_id, group in long_df.groupby("natural_protein_id", sort=True):
        per_experiment = group.groupby("experiment_id")
        q = {
            "natural_protein_id": protein_id,
            "protein_name": group["protein_name"].iloc[0],
            "n_rows": int(len(group)),
            "n_unique_8mers": int(group["dna_sequence"].nunique()),
            "n_rc_classes": int(group["canonical_rc"].nunique()),
            "n_missing_8mers": int(EXPECTED_8MER_RC_CLASSES - group["canonical_rc"].nunique()),
            "n_duplicate_rows": int(len(group) - group["dna_sequence"].nunique()),
            "n_experiments": int(per_experiment.ngroups),
            "n_replicate_groups": int(group["replicate_id"].nunique()),
            "dna_qc_pass": bool(group["dna_valid"].all()),
            "rc_qc_pass": bool(group["rc_source_matches"].all()),
            "quality_level": "usable" if group["dna_valid"].all() else "exclude",
        }
        qc_rows.append(q)
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(TABLES / "natural_pbm_qc_summary.csv", index=False)
    return qc


def replicate_qc(long_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein_id, group in long_df.groupby("natural_protein_id", sort=True):
        if group["experiment_id"].nunique() < 2:
            rows.append(
                {
                    "natural_protein_id": protein_id,
                    "protein_name": group["protein_name"].iloc[0],
                    "score_type": group["experimental_score_type"].iloc[0],
                    "n_experiments": int(group["experiment_id"].nunique()),
                    "pearson": np.nan,
                    "spearman": np.nan,
                    "replicate_status": "single_measurement",
                }
            )
            continue
        # Use first two experiments when a publication exposes replicate-like comparisons.
        exp_ids = list(group["experiment_id"].dropna().unique())[:2]
        a = group[group["experiment_id"] == exp_ids[0]][["canonical_rc", "experimental_score"]].rename(
            columns={"experimental_score": "score_a"}
        )
        b = group[group["experiment_id"] == exp_ids[1]][["canonical_rc", "experimental_score"]].rename(
            columns={"experimental_score": "score_b"}
        )
        merged = a.merge(b, on="canonical_rc", how="inner")
        if merged.empty:
            pearson = spearman = np.nan
        else:
            pearson = float(merged["score_a"].corr(merged["score_b"], method="pearson"))
            spearman = float(merged["score_a"].corr(merged["score_b"], method="spearman"))
        rows.append(
            {
                "natural_protein_id": protein_id,
                "protein_name": group["protein_name"].iloc[0],
                "score_type": group["experimental_score_type"].iloc[0],
                "n_experiments": int(group["experiment_id"].nunique()),
                "pearson": pearson,
                "spearman": spearman,
                "replicate_status": "replicated" if group["experiment_id"].nunique() > 1 else "single_measurement",
            }
        )
    rep = pd.DataFrame(rows)
    rep.to_csv(TABLES / "natural_pbm_replicate_qc.csv", index=False)
    return rep


def natural_benchmark(long_df: pd.DataFrame, meta: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    consensus = consensus_by_experiment(long_df)
    merged = consensus.merge(meta, on=["natural_protein_id", "protein_name"], how="left", suffixes=("", "_meta"))
    qc_subset = qc[["natural_protein_id", "quality_level"]]
    merged = merged.merge(qc_subset, on="natural_protein_id", how="left")
    merged = merged.rename(
        columns={
            "natural_protein_id": "protein_id",
            "dna_sequence": "dna_sequence",
            "dna_length": "dna_length",
            "experimental_score": "experimental_score",
            "experimental_score_type": "experimental_score_type",
            "publication": "source",
        }
    )
    merged["experimental_percentile"] = merged.groupby("protein_id")["experimental_score"].rank(pct=True, ascending=True)
    merged["experimental_rank"] = merged.groupby("protein_id")["experimental_score"].rank(method="average", ascending=False)
    pub_col = next((c for c in ["publication", "publication_x", "publication_y", "archive_publication"] if c in merged.columns), None)
    if pub_col is not None:
        merged["source"] = merged["source"].fillna(merged[pub_col])
    else:
        merged["source"] = merged["source"].fillna("UniPROBE")
    columns = [
        "protein_id",
        "protein_name",
        "protein_family",
        "species",
        "dna_sequence",
        "dna_length",
        "canonical_rc",
        "experimental_score",
        "experimental_percentile",
        "experimental_rank",
        "experiment_id",
        "assay_type",
        "score_type",
        "quality_level",
        "source",
        "sequence_type",
        "sequence_match_to_assay",
        "replicate_group",
        "n_rows_collapsed",
        "experiment_ids",
        "replicate_ids",
    ]
    out = merged.loc[:, [col for col in columns if col in merged.columns]].copy()
    out.to_parquet(PROCESSED / "natural_pbm_benchmark_v0_4_1.parquet", index=False)
    out.to_csv(TABLES / "natural_pbm_benchmark_v0_4_1.csv", index=False)
    return out


def main() -> None:
    long_df = build_metadata_frame()
    meta = protein_metadata(long_df)
    qc = quality_tables(long_df)
    rep = replicate_qc(long_df)
    bench = natural_benchmark(long_df, meta, qc)
    print(f"long rows={len(long_df)} proteins={meta.shape[0]} benchmark_rows={len(bench)}")
    print(rep.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
