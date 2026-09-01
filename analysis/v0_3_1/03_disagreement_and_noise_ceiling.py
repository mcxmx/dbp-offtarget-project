from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
PROCESSED_DIR = ROOT / "data" / "processed" / "v0_3_1"
TABLES_DIR = ensure_dir(ROOT / "results" / "v0_3_1" / "tables")
FIGURES_DIR = ensure_dir(ROOT / "results" / "v0_3_1" / "figures")
DOCS_DIR = ensure_dir(ROOT / "docs" / "v0_3_1")
V03_TABLES_DIR = ROOT / "results" / "v0_3" / "tables"
HIGH_SCORE_QUANTILE = 0.95
SIMILARITY_REFERENCE_QUANTILE = 0.50
TOP_EXAMPLES_PER_PROTEIN = 20


def add_percentiles(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein_id, group in scored.groupby("protein_id"):
        out = group.copy()
        out["experimental_percentile"] = out["experimental_escore_consensus"].rank(method="average", pct=True)
        out["experimental_rank"] = out["experimental_escore_consensus"].rank(method="average", ascending=False)
        rows.append(out)
    return pd.concat(rows, ignore_index=True)


def disagreement_tables(scored: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = []
    examples = []
    metric = "hamming_similarity_to_paper_motif_rc_aware"
    for protein_id, group in scored.groupby("protein_id"):
        score_threshold = float(group["experimental_escore_consensus"].quantile(HIGH_SCORE_QUANTILE))
        similarity_threshold = float(group[metric].quantile(SIMILARITY_REFERENCE_QUANTILE))
        candidates = group[
            (group["experimental_escore_consensus"] >= score_threshold)
            & (group[metric] <= similarity_threshold)
        ].copy()
        threshold_definition = (
            f"per-protein E-score >= {HIGH_SCORE_QUANTILE:.2f} quantile "
            f"and RC-aware motif hamming similarity <= {SIMILARITY_REFERENCE_QUANTILE:.2f} quantile"
        )
        counts.append(
            {
                "protein_id": protein_id,
                "n_total_sequences": int(len(group)),
                "n_disagreement": int(len(candidates)),
                "fraction_disagreement": float(len(candidates) / len(group)),
                "experimental_score_threshold": score_threshold,
                "sequence_similarity_threshold": similarity_threshold,
                "threshold_definition": threshold_definition,
                "notes": "Sequence-vs-experiment disagreement candidates; not model failures.",
            }
        )
        top = candidates.sort_values("experimental_escore_consensus", ascending=False).head(TOP_EXAMPLES_PER_PROTEIN).copy()
        top["example_limit_per_protein"] = TOP_EXAMPLES_PER_PROTEIN
        top["threshold_definition"] = threshold_definition
        top["notes"] = "Examples only; do not interpret row count in this file as total discoveries."
        examples.append(top)
    return pd.DataFrame(counts), pd.concat(examples, ignore_index=True) if examples else pd.DataFrame()


def noise_ceiling() -> pd.DataFrame:
    replicate = pd.read_csv(V03_TABLES_DIR / "replicate_qc.csv")
    rows = []
    for _, row in replicate[replicate["qc_status"] == "replicate_pair"].iterrows():
        rows.append(
            {
                "protein_id": row["protein_id"],
                "score_type": row["score_type"],
                "pearson_correlation": row["pearson_correlation"],
                "spearman_correlation": row["spearman_correlation"],
                "n_aligned_7mers": row["n_aligned_7mers"],
                "interpretation": "empirical replicate agreement / assay reproducibility reference; not a strict mathematical ceiling",
            }
        )
    return pd.DataFrame(rows)


def plot_noise(noise: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    sub = noise[noise["score_type"] == "e_score"].copy()
    long = sub.melt(
        id_vars=["protein_id", "score_type"],
        value_vars=["pearson_correlation", "spearman_correlation"],
        var_name="correlation_type",
        value_name="correlation",
    )
    fig, ax = plt.subplots(figsize=(8, 4.6))
    sns.barplot(data=long, x="protein_id", y="correlation", hue="correlation_type", ax=ax)
    ax.set_ylim(0, 1)
    ax.set_xlabel("designed DBP")
    ax.set_ylabel("replicate correlation")
    ax.set_title("uPBM E-score replicate agreement reference")
    ax.legend(title="")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig_noise_ceiling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_disagreement_definition() -> None:
    text = f"""# Sequence-vs-Experiment Disagreement Definition

Audit date: 2026-09-01

The v0.3 file `sequence_similarity_failure_candidates.csv` contained 140 rows because the analysis selected at most 20 example rows per protein. That number is an example-table size, not the total number of discoveries.

## v0.3.1 Terminology

Use: `sequence-vs-experiment disagreement cases`.

Do not use: `model failures`.

Reason: no protein-conditioned model has been evaluated yet.

## v0.3.1 Criterion

For each protein separately, using 8192 reverse-complement equivalence classes:

- Experimental high-score condition: processed uPBM E-score consensus >= the per-protein {HIGH_SCORE_QUANTILE:.0%} quantile.
- Sequence-similarity condition: RC-aware hamming similarity to the paper motif <= the per-protein {SIMILARITY_REFERENCE_QUANTILE:.0%} quantile.

This is a percentile-based descriptive rule. It is not an assay-defined binding threshold and should not be interpreted as an absolute off-target risk threshold.

## Outputs

- `results/v0_3_1/tables/all_disagreement_candidate_counts.csv`: total counts by protein.
- `results/v0_3_1/tables/top_disagreement_examples.csv`: at most {TOP_EXAMPLES_PER_PROTEIN} examples per protein, for inspection only.
"""
    (DOCS_DIR / "DISAGREEMENT_DEFINITION.md").write_text(text, encoding="utf-8")


def main() -> None:
    scored = pd.read_parquet(PROCESSED_DIR / "designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1.parquet")
    scored = add_percentiles(scored)
    counts, examples = disagreement_tables(scored)
    counts.to_csv(TABLES_DIR / "all_disagreement_candidate_counts.csv", index=False)
    examples.to_csv(TABLES_DIR / "top_disagreement_examples.csv", index=False)
    noise = noise_ceiling()
    noise.to_csv(TABLES_DIR / "experimental_noise_ceiling.csv", index=False)
    plot_noise(noise)
    write_disagreement_definition()
    print(counts.to_string(index=False))
    print(noise[noise["score_type"] == "e_score"].to_string(index=False))
    print(f"top examples written: {len(examples)}")


if __name__ == "__main__":
    main()
