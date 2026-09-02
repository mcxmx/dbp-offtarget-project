from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
RESULTS = ensure_dir(ROOT / "results" / "v0_4_1")
TABLES = ensure_dir(RESULTS / "tables")
FIGURES = ensure_dir(RESULTS / "figures")
DOCS = ensure_dir(ROOT / "docs" / "v0_4_1")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4_1")
TODAY = date.today().isoformat()


def unique_semicolon_count(series: pd.Series) -> int:
    items = set()
    for value in series.dropna().astype(str):
        for part in value.split(";"):
            if part:
                items.add(part)
    return len(items)


def read_metric(macro: pd.DataFrame, baseline: str, metric: str, dataset: str | None = None) -> tuple[float, int]:
    sub = macro[(macro["baseline"] == baseline) & (macro["metric"] == metric)].copy()
    if dataset is not None and "dataset" in sub.columns:
        sub = sub[sub["dataset"] == dataset]
    if sub.empty:
        return np.nan, 0
    row = sub.iloc[0]
    n_col = "n_proteins" if "n_proteins" in row.index else "n_proteins_with_metric"
    return float(row["median"]) if pd.notna(row["median"]) else np.nan, int(row.get(n_col, 0))


def main() -> None:
    natural = pd.read_parquet(ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet")
    qc = pd.read_csv(TABLES / "natural_pbm_qc_summary.csv")
    rep = pd.read_csv(TABLES / "natural_pbm_replicate_qc.csv")
    splits = pd.read_csv(METADATA / "natural_pbm_splits.csv")
    clusters = pd.read_csv(METADATA / "natural_protein_clusters.csv")
    spc_perf = pd.read_csv(TABLES / "simple_pc_performance.csv")
    spc_macro = pd.read_csv(TABLES / "simple_pc_performance_macro.csv")
    seq_macro = pd.read_csv(ROOT / "results" / "v0_4" / "tables" / "baseline_performance_macro.csv")
    dis = pd.read_csv(TABLES / "simple_pc_disagreement_resolution.csv")
    deeppbs_perf = pd.read_csv(TABLES / "deeppbs_performance.csv")

    seq_median, seq_n = read_metric(seq_macro, "sequence_kmer3", "spearman")
    nampnn_median, nampnn_n = read_metric(seq_macro, "NA-MPNN_structural_ppm", "spearman")
    spc_nat, spc_nat_n = read_metric(spc_macro, "SimpleProteinConditionalBaseline_composition_ridge", "spearman", "natural_test")
    spc_des, spc_des_n = read_metric(spc_macro, "SimpleProteinConditionalBaseline_composition_ridge", "spearman", "designed_external")
    replicate_ref = 0.591
    resolved_total = int(dis["n_resolved"].sum())
    disagreement_total = int(dis["n_v0_3_1_disagreement_candidates"].sum())
    resolution_rate = resolved_total / disagreement_total
    split_summary = splits.groupby("split")["protein_id"].nunique().to_dict()

    summary_rows = [
        {
            "baseline": "sequence_kmer3_RC_aware",
            "training_source": "none",
            "evaluation_dataset": "GSE237017 designed uPBM",
            "coverage_proteins": seq_n,
            "macro_median_spearman": seq_median,
            "notes": "Best v0.3.1/v0.4 sequence-only baseline; not protein-conditioned.",
        },
        {
            "baseline": "SimpleProteinConditionalBaseline_composition_ridge",
            "training_source": "UniPROBE natural PBM train split",
            "evaluation_dataset": "UniPROBE natural PBM natural_test split",
            "coverage_proteins": spc_nat_n,
            "macro_median_spearman": spc_nat,
            "notes": "Lightweight protein-conditioned baseline using reference full-length UniProt sequences; not proposed method.",
        },
        {
            "baseline": "SimpleProteinConditionalBaseline_composition_ridge",
            "training_source": "UniPROBE natural PBM train split",
            "evaluation_dataset": "GSE237017 designed uPBM external",
            "coverage_proteins": spc_des_n,
            "macro_median_spearman": spc_des,
            "notes": "Natural-to-designed external test with assay/k-mer-length caveats.",
        },
        {
            "baseline": "NA-MPNN_structural_ppm",
            "training_source": "official checkpoint / diagnostic v0.4",
            "evaluation_dataset": "GSE237017 designed uPBM",
            "coverage_proteins": nampnn_n,
            "macro_median_spearman": nampnn_median,
            "notes": "Diagnostic only; 2/7 proteins covered and DBP48 has overlap risk.",
        },
        {
            "baseline": "DeepPBS",
            "training_source": "official checkpoint intended",
            "evaluation_dataset": "GSE237017 designed uPBM",
            "coverage_proteins": int((deeppbs_perf["n_rc_classes"] > 0).sum()),
            "macro_median_spearman": np.nan,
            "notes": "Not run in v0.4.1 because this host lacks Docker/WSL runtime for official Linux workflow.",
        },
        {
            "baseline": "designed_uPBM_replicate_reference",
            "training_source": "not a model",
            "evaluation_dataset": "GSE237017 designed uPBM replicate agreement",
            "coverage_proteins": 7,
            "macro_median_spearman": replicate_ref,
            "notes": "Empirical assay reproducibility reference from v0.3.1, not a strict theoretical upper bound.",
        },
    ]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TABLES / "v0_4_1_baseline_summary.csv", index=False)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(7.5, 4))
    plot_df = summary[summary["baseline"] != "DeepPBS"].copy()
    ax.bar(plot_df["baseline"] + "\n" + plot_df["evaluation_dataset"], plot_df["macro_median_spearman"], color="#3A6EA5")
    ax.axhline(replicate_ref, color="black", linestyle="--", linewidth=1, label="replicate reference")
    ax.set_ylabel("Macro median Spearman")
    ax.set_title("Available Baselines vs Experimental Replicate Reference")
    ax.tick_params(axis="x", labelrotation=35, labelsize=7)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_v0_4_1_6_all_baselines_vs_replicate_reference.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 3.8))
    dis_plot = dis[["protein_id", "n_resolved", "n_v0_3_1_disagreement_candidates"]].copy()
    dis_plot["n_unresolved"] = dis_plot["n_v0_3_1_disagreement_candidates"] - dis_plot["n_resolved"]
    ax.bar(dis_plot["protein_id"], dis_plot["n_resolved"], label="resolved by SimplePC top-10%", color="#3A6EA5")
    ax.bar(dis_plot["protein_id"], dis_plot["n_unresolved"], bottom=dis_plot["n_resolved"], label="unresolved", color="#C95D63")
    ax.set_ylabel("Sequence-vs-experiment disagreement candidates")
    ax.set_title("SimplePC Resolution of v0.3.1 Disagreements")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_v0_4_1_7_failure_case_resolution.png", dpi=300)
    plt.close(fig)

    (DOCS / "NATURAL_TO_DESIGNED_EVALUATION_PLAN.md").write_text(
        f"""# v0.4.1 Natural-to-Designed Evaluation Plan

Date: {TODAY}

## Train/Test Design

- Train: UniPROBE natural PBM train proteins from `metadata/v0_4_1/natural_pbm_splits.csv`.
- Validation: UniPROBE natural PBM validation proteins; no designed DBP rows used for hyperparameter choice.
- Test A: held-out natural PBM proteins.
- Test B: GSE237017 designed DBPs, kept external.

## Confounders

- Natural UniPROBE uses processed contiguous 8-mer E-scores.
- Designed GSE237017 uses processed uPBM 7-mer E-scores.
- Protein sequences in v0.4.1 are full-length UniProt references, not confirmed assay constructs.
- A natural-to-designed performance drop may reflect assay/platform/k-mer processing shift as well as protein-distribution shift.

## Split Rules

- Protein split is cluster-aware at the 40% proxy identity level.
- DNA rows are grouped by reverse-complement canonical class; no oriented RC pair is split into separate evaluation units.
- Designed DBP sequences are excluded from natural training.
""",
        encoding="utf-8",
    )

    gate = "WAIT FOR STRONGER BASELINE"
    gate_text = f"""# v0.4.1 Final Model Development Gate

Date: {TODAY}

## Decision

{gate}

## Evidence

1. Natural PBM benchmark: {natural['protein_id'].nunique():,} proteins and {len(natural):,} protein-RC-class 8-mer units are now available for simple training.
2. Natural held-out SimplePC macro median Spearman: {spc_nat:.3f}.
3. Designed external SimplePC macro median Spearman: {spc_des:.3f}.
4. Best prior sequence-only designed baseline macro median Spearman: {seq_median:.3f}.
5. Designed uPBM empirical replicate Spearman reference: about {replicate_ref:.3f}.
6. DeepPBS was not fairly run in v0.4.1 because this host lacks Docker/WSL runtime for the official Linux preprocessing workflow.
7. NA-MPNN remains diagnostic only: {nampnn_n}/7 designed proteins covered; DBP48 has known overlap risk.
8. SimplePC resolves {resolved_total:,}/{disagreement_total:,} v0.3.1 sequence-vs-experiment disagreement candidates ({resolution_rate:.1%}), leaving most unresolved.

## Interpretation

SimplePC improves over the sequence-only baseline on designed DBPs, but it remains well below the empirical replicate reference and fails on DBP6/DBP48. Because DeepPBS has not yet been fairly reproduced, the project should not start the final proposed model as the next step. The immediate next baseline task is a Linux/Docker DeepPBS run on all structurally evaluable designed DBPs, followed by the same per-protein ranking evaluation.
"""
    (RESULTS / "FINAL_MODEL_DEVELOPMENT_GATE.md").write_text(gate_text, encoding="utf-8")

    (DOCS / "PROPOSED_MODEL_REQUIREMENTS_V2.md").write_text(
        f"""# Proposed Model Requirements V2

Date: {TODAY}

These are conditional requirements, not an implemented architecture.

## Requirements Driven by v0.4.1

1. Support variable-length DNA k-mers at minimum 7-mer and 8-mer inputs without cropping/padding that changes the biological task.
2. Use protein information beyond amino-acid composition, because the SimplePC composition ridge still leaves a large gap to replicate agreement.
3. Preserve per-protein ranking evaluation and RC canonicalization as first-class constraints.
4. Explicitly model natural-to-designed generalization without conflating assay shift with biological OOD.
5. Include structure-aware baseline comparison once DeepPBS is runnable in Linux/Docker.

The final model should not be started until DeepPBS is fairly reproduced or formally ruled out as non-evaluable for the designed benchmark.
""",
        encoding="utf-8",
    )

    report = f"""# v0.4.1 Summary Report

Date: {TODAY}

## Natural PBM

- Final natural proteins: {natural['protein_id'].nunique():,}
- Coarse protein families: {natural[['protein_id', 'protein_family']].drop_duplicates()['protein_family'].nunique():,}
- Species: {natural[['protein_id', 'species']].drop_duplicates()['species'].nunique():,}
- Experiments/profile groups: {unique_semicolon_count(natural['experiment_id']):,}
- Protein-DNA units: {len(natural):,}
- Assay: UniPROBE universal PBM, processed contiguous 8-mer E-score.
- K-mer length: 8 bp RC classes.
- Protein sequence completeness in final benchmark: 100%; construct sequence completeness remains 0%.
- Replicated protein/construct groups in downloaded metadata: {int(rep['spearman'].notna().sum()):,}; median Spearman {rep['spearman'].dropna().median():.3f} where available.
- 40% cluster split: train {split_summary.get('train', 0)}, validation {split_summary.get('validation', 0)}, natural_test {split_summary.get('natural_test', 0)} proteins.

## Baselines

- Sequence-only designed median Spearman: {seq_median:.3f}.
- SimplePC natural held-out median Spearman: {spc_nat:.3f}.
- SimplePC designed external median Spearman: {spc_des:.3f}.
- NA-MPNN diagnostic designed median Spearman: {nampnn_median:.3f} over {nampnn_n} proteins.
- DeepPBS: not run; Docker/WSL unavailable on this host.

## Failure Summary

SimplePC improves designed external median Spearman over the best sequence-only metric by {spc_des - seq_median:.3f}, but remains {replicate_ref - spc_des:.3f} below the empirical replicate reference. It resolves {resolved_total:,}/{disagreement_total:,} pre-registered sequence-vs-experiment disagreement candidates. DBP6 and DBP48 remain the clearest designed-protein failures for this baseline.
"""
    (RESULTS / "V0_4_1_SUMMARY_REPORT.md").write_text(report, encoding="utf-8")
    print(report)
    print(gate_text)


if __name__ == "__main__":
    main()
