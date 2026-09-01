from __future__ import annotations

import itertools
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, load_yaml, project_root, reverse_complement


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
EXPECTED_7MERS = int(CONFIG["benchmark_v0_3"]["expected_7mer_count"])
INTERIM_DIR = ensure_dir(ROOT / "data" / "interim" / "gse237017")
METADATA_DIR = ensure_dir(ROOT / "metadata" / "v0_3")
RESULTS_DIR = ensure_dir(ROOT / "results" / "v0_3")
FIGURES_DIR = ensure_dir(RESULTS_DIR / "figures")
TABLES_DIR = ensure_dir(RESULTS_DIR / "tables")

SCORE_COLUMNS = ["e_score", "median_intensity", "z_score"]


def spearman(x: pd.Series, y: pd.Series) -> float:
    return x.rank(method="average").corr(y.rank(method="average"), method="pearson")


def replicate_qc(df: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein_id, sample_group in samples.groupby("protein_id"):
        sample_ids = sample_group.sort_values(["replicate", "gsm_id"])["gsm_id"].tolist()
        if len(sample_ids) < 2:
            for score_col in SCORE_COLUMNS:
                rows.append(
                    {
                        "protein_id": protein_id,
                        "gsm_a": sample_ids[0],
                        "replicate_a": sample_group.iloc[0]["replicate"],
                        "gsm_b": "",
                        "replicate_b": "",
                        "score_type": score_col,
                        "n_aligned_7mers": 0,
                        "pearson_correlation": pd.NA,
                        "spearman_correlation": pd.NA,
                        "qc_status": "single_replicate_only",
                        "notes": "No replicate correlation calculated because GEO metadata contains one sample for this protein.",
                    }
                )
            continue
        for gsm_a, gsm_b in itertools.combinations(sample_ids, 2):
            a = df[df["gsm_id"] == gsm_a][["dna_7mer", *SCORE_COLUMNS]].copy()
            b = df[df["gsm_id"] == gsm_b][["dna_7mer", *SCORE_COLUMNS]].copy()
            merged = a.merge(b, on="dna_7mer", suffixes=("_a", "_b"), how="inner")
            replicate_a = samples.loc[samples["gsm_id"] == gsm_a, "replicate"].iloc[0]
            replicate_b = samples.loc[samples["gsm_id"] == gsm_b, "replicate"].iloc[0]
            for score_col in SCORE_COLUMNS:
                left = merged[f"{score_col}_a"]
                right = merged[f"{score_col}_b"]
                valid = left.notna() & right.notna()
                if valid.sum() == 0:
                    pearson = pd.NA
                    spear = pd.NA
                else:
                    pearson = left[valid].corr(right[valid], method="pearson")
                    spear = spearman(left[valid], right[valid])
                rows.append(
                    {
                        "protein_id": protein_id,
                        "gsm_a": gsm_a,
                        "replicate_a": replicate_a,
                        "gsm_b": gsm_b,
                        "replicate_b": replicate_b,
                        "score_type": score_col,
                        "n_aligned_7mers": int(valid.sum()),
                        "pearson_correlation": pearson,
                        "spearman_correlation": spear,
                        "qc_status": "replicate_pair",
                        "notes": "Replicates aligned by expanded dna_7mer.",
                    }
                )
    return pd.DataFrame(rows).sort_values(["protein_id", "score_type", "gsm_a", "gsm_b"])


def reverse_complement_qc(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (protein_id, gsm_id), group in df.groupby(["protein_id", "gsm_id"]):
        score_map = group.set_index("dna_7mer")[SCORE_COLUMNS]
        sample_rows = {"protein_id": protein_id, "gsm_id": gsm_id, "n_7mers": int(len(score_map))}
        for score_col in SCORE_COLUMNS:
            diffs = []
            for seq, score in score_map[score_col].items():
                rc = reverse_complement(seq)
                if rc in score_map.index:
                    diffs.append(abs(score - score_map.loc[rc, score_col]))
            series = pd.Series(diffs, dtype="float64")
            sample_rows[f"{score_col}_n_pairs_checked"] = int(len(series))
            sample_rows[f"{score_col}_max_abs_rc_difference"] = float(series.max()) if not series.empty else pd.NA
            sample_rows[f"{score_col}_mean_abs_rc_difference"] = float(series.mean()) if not series.empty else pd.NA
        sample_rows["notes"] = "Processed table contains paired 7-mer and reverse-complement columns; parser expands both columns and checks score symmetry."
        rows.append(sample_rows)
    return pd.DataFrame(rows).sort_values(["protein_id", "gsm_id"])


def plot_sample_overview(samples: pd.DataFrame, coverage: pd.DataFrame) -> None:
    overview = (
        coverage.groupby("protein_id")
        .agg(
            n_samples=("gsm_id", "nunique"),
            min_unique_7mers=("n_unique_7mers", "min"),
            max_unique_7mers=("n_unique_7mers", "max"),
        )
        .reset_index()
    )
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    sns.barplot(data=overview, x="protein_id", y="n_samples", ax=axes[0], color="#4c72b0")
    axes[0].set_xlabel("designed DBP")
    axes[0].set_ylabel("number of GEO samples")
    axes[0].set_title("uPBM sample count")
    sns.barplot(data=overview, x="protein_id", y="min_unique_7mers", ax=axes[1], color="#55a868")
    axes[1].axhline(EXPECTED_7MERS, linestyle="--", color="black", linewidth=1, label="expected 4^7")
    axes[1].set_xlabel("designed DBP")
    axes[1].set_ylabel("unique 7-mers per sample")
    axes[1].set_title("7-mer coverage")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_v0_3_1_dataset_sample_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_replicate_correlations(rep_qc: pd.DataFrame) -> None:
    plot_df = rep_qc[rep_qc["qc_status"] == "replicate_pair"].copy()
    plot_df["pearson_correlation"] = pd.to_numeric(plot_df["pearson_correlation"], errors="coerce")
    sns.set_theme(style="whitegrid", context="paper")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    sns.barplot(data=plot_df, x="protein_id", y="pearson_correlation", hue="score_type", ax=ax)
    ax.set_xlabel("designed DBP")
    ax.set_ylabel("replicate Pearson correlation")
    ax.set_title("uPBM replicate consistency")
    ax.set_ylim(0, 1)
    ax.legend(title="score")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_v0_3_2_replicate_correlations.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "replicate_correlations.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    df = pd.read_parquet(INTERIM_DIR / "upbm_7mers_long.parquet")
    samples = pd.read_csv(METADATA_DIR / "gse237017_samples.csv", dtype={"replicate": str})
    coverage = pd.read_csv(TABLES_DIR / "sample_coverage_qc.csv", dtype={"replicate": str})
    rep_qc = replicate_qc(df, samples)
    rc_qc = reverse_complement_qc(df)
    rep_qc.to_csv(TABLES_DIR / "replicate_qc.csv", index=False)
    rc_qc.to_csv(TABLES_DIR / "reverse_complement_qc.csv", index=False)
    plot_sample_overview(samples, coverage)
    plot_replicate_correlations(rep_qc)
    valid = rep_qc[rep_qc["qc_status"] == "replicate_pair"].copy()
    print(valid.groupby("score_type")[["pearson_correlation", "spearman_correlation"]].agg(["min", "median", "max"]).to_string())
    print("single replicate only:", sorted(rep_qc.loc[rep_qc["qc_status"] == "single_replicate_only", "protein_id"].unique()))


if __name__ == "__main__":
    main()
