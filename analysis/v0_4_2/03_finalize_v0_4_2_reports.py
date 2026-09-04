from __future__ import annotations

from datetime import date
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
TODAY = date.today().isoformat()
V042_RESULTS = ensure_dir(ROOT / "results" / "v0_4_2")
V042_TABLES = ensure_dir(V042_RESULTS / "tables")
V042_FIGURES = ensure_dir(V042_RESULTS / "figures")
V042_DOCS = ensure_dir(ROOT / "docs" / "v0_4_2")


def main() -> None:
    natural_cov = pd.read_csv(V042_TABLES / "natural_construct_coverage.csv")
    frozen_macro = pd.read_csv(V042_TABLES / "frozen_plm_performance_macro.csv")
    simple_pc = pd.read_csv(ROOT / "results" / "v0_4_1" / "tables" / "simple_pc_performance_macro.csv")
    seq_macro = pd.read_csv(ROOT / "results" / "v0_4" / "tables" / "baseline_performance_macro.csv")
    replicate = pd.read_csv(ROOT / "results" / "v0_3_1" / "tables" / "experimental_noise_ceiling.csv")
    difficulty = pd.read_csv(V042_TABLES / "designed_difficulty_factors.csv")
    resolution = pd.read_csv(V042_TABLES / "disagreement_resolution_v0_4_2.csv")
    disagreement_counts = pd.read_csv(
        ROOT / "results" / "v0_3_1" / "tables" / "all_disagreement_candidate_counts.csv"
    )
    common_hard = pd.read_parquet(V042_RESULTS.parent.parent / "data" / "processed" / "v0_4_2" / "common_hard_specificity_cases.parquet")
    disagreement_total = int(disagreement_counts["n_disagreement"].sum())
    common_high = common_hard[common_hard["failure_type"].eq("high_experiment_low_all_core")]
    resolution_core = (
        resolution[
            resolution["method"].isin(
                [
                    "sequence_kmer3",
                    "SimpleProteinConditionalBaseline",
                    "FrozenPLMProteinConditionalBaseline",
                ]
            )
        ]
        .groupby("method", as_index=False)[["n_evaluable", "n_resolved"]]
        .sum()
    )
    hardest = difficulty.loc[difficulty["frozen_plm_spearman"].idxmin()]

    summary_rows = [
        {
            "method": "sequence-only best",
            "protein_representation": "none",
            "structure_required": False,
            "training_data": "v0.2/v0.1 proxy benchmark",
            "n_designed_proteins_covered": 7,
            "designed_macro_median_spearman": float(seq_macro[seq_macro["baseline"].eq("sequence_kmer3")]["median"].iloc[0]),
            "natural_macro_median_spearman": None,
            "gap_to_replicate": None,
            "overlap_caveat": "sequence-only proxy baseline",
            "notes": "best prior designed sequence-only baseline",
        },
        {
            "method": "SimpleProteinConditionalBaseline",
            "protein_representation": "simple composition features",
            "structure_required": False,
            "training_data": "natural PBM train split",
            "n_designed_proteins_covered": 7,
            "designed_macro_median_spearman": float(simple_pc[simple_pc["dataset"].eq("designed_external")]["median"].iloc[0]),
            "natural_macro_median_spearman": float(simple_pc[simple_pc["dataset"].eq("natural_test")]["median"].iloc[0]),
            "gap_to_replicate": float(replicate[replicate["score_type"].eq("e_score")]["spearman_correlation"].median()) - float(simple_pc[simple_pc["dataset"].eq("designed_external")]["median"].iloc[0]),
            "overlap_caveat": "no structure requirement; designed test external",
            "notes": "previous baseline",
        },
        {
            "method": "FrozenPLMProteinConditionalBaseline",
            "protein_representation": "frozen ESM-2 mean-pooled embeddings",
            "structure_required": False,
            "training_data": "natural PBM train split",
            "n_designed_proteins_covered": 7,
            "designed_macro_median_spearman": float(frozen_macro[frozen_macro["dataset"].eq("designed_external")]["median"].iloc[0]),
            "natural_macro_median_spearman": float(frozen_macro[frozen_macro["dataset"].eq("natural_test")]["median"].iloc[0]),
            "gap_to_replicate": float(replicate[replicate["score_type"].eq("e_score")]["spearman_correlation"].median()) - float(frozen_macro[frozen_macro["dataset"].eq("designed_external")]["median"].iloc[0]),
            "overlap_caveat": "no structure requirement; designed test external",
            "notes": "frozen protein-LM baseline",
        },
        {
            "method": "DeepPBS",
            "protein_representation": "structure-aware",
            "structure_required": True,
            "training_data": "official model provenance only",
            "n_designed_proteins_covered": 0,
            "designed_macro_median_spearman": None,
            "natural_macro_median_spearman": None,
            "gap_to_replicate": None,
            "overlap_caveat": "official Linux example not run on host",
            "notes": "environment-limited in this repository host",
        },
        {
            "method": "NA-MPNN diagnostic",
            "protein_representation": "structure-aware",
            "structure_required": True,
            "training_data": "official diagnostic run",
            "n_designed_proteins_covered": 2,
            "designed_macro_median_spearman": float(pd.read_csv(ROOT / "results" / "v0_4" / "tables" / "baseline_performance_macro.csv").query("baseline == 'NA-MPNN_structural_ppm' and metric == 'spearman'")["median"].iloc[0]),
            "natural_macro_median_spearman": None,
            "gap_to_replicate": None,
            "overlap_caveat": "diagnostic only; DBP48 overlap risk",
            "notes": "not a full strong-baseline result",
        },
        {
            "method": "Replicate reference",
            "protein_representation": "n/a",
            "structure_required": False,
            "training_data": "experimental agreement reference",
            "n_designed_proteins_covered": 7,
            "designed_macro_median_spearman": float(replicate[replicate["score_type"].eq("e_score")]["spearman_correlation"].median()),
            "natural_macro_median_spearman": None,
            "gap_to_replicate": 0.0,
            "overlap_caveat": "empirical assay reproducibility reference",
            "notes": "not a theoretical upper bound",
        },
    ]
    pd.DataFrame(summary_rows).to_csv(V042_TABLES / "final_strong_baseline_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    bars = pd.DataFrame(summary_rows).set_index("method")
    plot_methods = [
        "sequence-only best",
        "SimpleProteinConditionalBaseline",
        "FrozenPLMProteinConditionalBaseline",
        "Replicate reference",
    ]
    values = [bars.loc[m, "designed_macro_median_spearman"] for m in plot_methods]
    colors = ["#4C78A8", "#54A24B", "#F58518", "#000000"]
    ax.bar(plot_methods, values, color=colors)
    ax.set_ylabel("Designed external macro median Spearman")
    ax.set_title("v0.4.2 Baseline Overview")
    ax.set_ylim(0, max(v for v in values if pd.notna(v)) * 1.15)
    ax.tick_params(axis="x", labelrotation=20)
    fig.tight_layout()
    fig.savefig(V042_FIGURES / "fig_v0_4_2_1_baseline_overview.png", dpi=300)
    plt.close(fig)

    gate_text = f"""# v0.4.2 Final Strong Baseline Gate

Date: {TODAY}

## Decision

WAIT - BENCHMARK STILL INCOMPLETE

## Evidence

1. Natural PBM construct audit recovered no assay-aligned construct sequences from the current local provenance for the 57 UniPROBE proteins.
2. FrozenPLMProteinConditionalBaseline improved the natural-test macro median Spearman to {float(frozen_macro[frozen_macro["dataset"].eq("natural_test")]["median"].iloc[0]):.3f}.
3. FrozenPLM designed-external macro median Spearman is {float(frozen_macro[frozen_macro["dataset"].eq("designed_external")]["median"].iloc[0]):.3f}.
4. The prior SimplePC designed-external macro median Spearman was {float(simple_pc[simple_pc["dataset"].eq("designed_external")]["median"].iloc[0]):.3f}.
5. The best prior designed sequence-only baseline remains {float(seq_macro[seq_macro["baseline"].eq("sequence_kmer3")]["median"].iloc[0]):.3f}.
6. Designed uPBM empirical replicate Spearman reference remains about {float(replicate[replicate["score_type"].eq("e_score")]["spearman_correlation"].median()):.3f}.
7. DeepPBS official Linux example was not run on this host because Docker/Podman/installed WSL are unavailable.
8. NA-MPNN remains diagnostic only; the prior v0.4 result covered 2/7 proteins.
9. The {disagreement_total:,} sequence-vs-experiment disagreement candidates from v0.3.1 are unchanged as the disagreement reference set.
10. Under the frozen resolution rule, the complete core methods resolve {int(resolution_core.loc[resolution_core["method"].eq("sequence_kmer3"), "n_resolved"].iloc[0]):,} (k-mer3), {int(resolution_core.loc[resolution_core["method"].eq("SimpleProteinConditionalBaseline"), "n_resolved"].iloc[0]):,} (SimplePC), and {int(resolution_core.loc[resolution_core["method"].eq("FrozenPLMProteinConditionalBaseline"), "n_resolved"].iloc[0]):,} (FrozenPLM) candidates.
11. The core common-hard set contains {len(common_high):,} high-experiment/low-all-core cases. This is a ranking discrepancy set, not a claim that every available model failed.
12. The lowest FrozenPLM designed Spearman is {hardest["frozen_plm_spearman"]:.3f} for {hardest["protein_id"]}; DBP6 and DBP48 remain difficult, but this descriptive result is not a biological mechanism conclusion.

## Interpretation

FrozenPLM is a real, frozen protein-conditioned baseline, but it did not improve designed-uPBM ranking relative to either the best sequence-only proxy or SimplePC. It remains far below empirical replicate agreement. Because the assay-aligned natural construct benchmark remains empty and DeepPBS is still not runnable on this host, the project should not yet transition to the final proposed model implementation from this repository state alone.
"""
    (V042_RESULTS / "FINAL_STRONG_BASELINE_GATE.md").write_text(gate_text, encoding="utf-8")

    report = f"""# v0.4.2 Validation Report

- Natural construct audit: {int(natural_cov["n_total_proteins"].iloc[0])} proteins audited; {int(natural_cov["n_unknown"].iloc[0])} unknown assay constructs.
- FrozenPLM designed macro median Spearman: {float(frozen_macro[frozen_macro["dataset"].eq("designed_external")]["median"].iloc[0]):.3f}.
- FrozenPLM natural macro median Spearman: {float(frozen_macro[frozen_macro["dataset"].eq("natural_test")]["median"].iloc[0]):.3f}.
- Best prior SimplePC designed macro median Spearman: {float(simple_pc[simple_pc["dataset"].eq("designed_external")]["median"].iloc[0]):.3f}.
- Best prior sequence-only designed macro median Spearman: {float(seq_macro[seq_macro["baseline"].eq("sequence_kmer3")]["median"].iloc[0]):.3f}.
- Empirical replicate e_score Spearman reference: {float(replicate[replicate["score_type"].eq("e_score")]["spearman_correlation"].median()):.3f}.
- {disagreement_total:,} disagreement candidates remain the reference disagreement set.
- Core common-hard high-experiment/low-all-core cases: {len(common_high):,}.
- Lowest FrozenPLM designed protein: {hardest["protein_id"]} ({hardest["frozen_plm_spearman"]:.3f} Spearman).

Overall decision: WAIT - BENCHMARK STILL INCOMPLETE.
"""
    (V042_RESULTS / "V0_4_2_VALIDATION_REPORT.md").write_text(report, encoding="utf-8")

    gap_report = f"""# v0.4.2 Baseline Gap Analysis

Date: {TODAY}

## Main ranking results

| Method | Designed macro-median Spearman | Natural-test macro-median Spearman | Coverage |
|---|---:|---:|---:|
| Best sequence-only k-mer3 proxy | {float(seq_macro[seq_macro["baseline"].eq("sequence_kmer3")]["median"].iloc[0]):.3f} | NA | 7/7 |
| SimpleProteinConditionalBaseline | {float(simple_pc[simple_pc["dataset"].eq("designed_external")]["median"].iloc[0]):.3f} | {float(simple_pc[simple_pc["dataset"].eq("natural_test")]["median"].iloc[0]):.3f} | 7/7 |
| FrozenPLMProteinConditionalBaseline | {float(frozen_macro[frozen_macro["dataset"].eq("designed_external")]["median"].iloc[0]):.3f} | {float(frozen_macro[frozen_macro["dataset"].eq("natural_test")]["median"].iloc[0]):.3f} | 7/7 / natural test |
| DeepPBS | NA | NA | 0/7 |
| NA-MPNN diagnostic | partial | NA | 2/7 |
| Empirical replicate reference | {float(replicate[replicate["score_type"].eq("e_score")]["spearman_correlation"].median()):.3f} | NA | 7/7 reference |

The replicate value is an empirical assay reproducibility reference, not a
strict theoretical upper bound. Scores are compared per protein before
macro-summary; raw scores are not pooled across proteins.

## Failure diagnostics

- The pre-registered v0.3.1 disagreement set contains {disagreement_total:,} RC-class candidates.
- Resolution counts under the v0.4.2 top-10% rule are recorded in
  `tables/disagreement_resolution_v0_4_2.csv`; missing predictions are
  `not_evaluable`, never zero-filled.
- The complete-core common-hard set contains {len(common_high):,} rows:
  experimental percentile >= 0.95 and sequence k-mer3, SimplePC, and FrozenPLM
  prediction percentiles <= 0.50.
- Common-hard rows are distributed across all seven DBPs:
  {", ".join(f"{protein}={int(count)}" for protein, count in common_high.groupby("protein_id").size().reindex(["DBP1", "DBP3", "DBP35", "DBP48", "DBP5", "DBP6", "DBP9"], fill_value=0).items())}.
- The lowest FrozenPLM per-protein result is {hardest["protein_id"]} at
  {hardest["frozen_plm_spearman"]:.3f} Spearman. DBP48 is also negative
  ({float(difficulty.loc[difficulty["protein_id"].eq("DBP48"), "frozen_plm_spearman"].iloc[0]):.3f}).
- Motif-distance stratification is descriptive and does not prove epistasis.
  The output is `tables/performance_by_motif_distance_v0_4_2.csv`.

## Interpretation

FrozenPLM provides a stronger protein representation than the composition
baseline on the natural test split, but its designed external median is below
both the SimplePC result and the best sequence-only k-mer3 result. This is a
real baseline result, not evidence for a proposed architecture. The current
gap analysis is insufficient for a final model-development GO because DeepPBS
was not executed in a Linux-compatible runtime and the assay-aligned natural
construct benchmark remains empty.
"""
    (V042_RESULTS / "BASELINE_GAP_ANALYSIS.md").write_text(gap_report, encoding="utf-8")

    progress = f"""# v0.4.2 Progress

Date: {TODAY}

## Completed

- Audited the natural UniPROBE PBM benchmark for assay construct provenance.
- Confirmed the local provenance does not recover assay-aligned construct sequences for the 57 natural proteins.
- Kept `FULL_LENGTH_REFERENCE` as a sensitivity benchmark only.
- Recorded DeepPBS official Linux runtime status and preserved Docker/WSL constraints in provenance.
- Trained and evaluated a frozen ESM-2 `esm2_t12_35M_UR50D` protein-conditioned baseline.
- Added difficulty-factor, training-space similarity, motif-distance, disagreement-resolution, and common-hard diagnostics.

## Current numbers

- Natural proteins audited: 57
- Assay-aligned construct benchmark rows: 0
- FrozenPLM natural-test macro median Spearman: {float(frozen_macro[frozen_macro["dataset"].eq("natural_test")]["median"].iloc[0]):.3f}
- FrozenPLM designed-external macro median Spearman: {float(frozen_macro[frozen_macro["dataset"].eq("designed_external")]["median"].iloc[0]):.3f}
- SimplePC designed-external macro median Spearman: {float(simple_pc[simple_pc["dataset"].eq("designed_external")]["median"].iloc[0]):.3f}
- Best prior sequence-only designed macro median Spearman: {float(seq_macro[seq_macro["baseline"].eq("sequence_kmer3")]["median"].iloc[0]):.3f}
- Empirical replicate Spearman reference: about {float(replicate[replicate["score_type"].eq("e_score")]["spearman_correlation"].median()):.3f}
- v0.3.1 disagreement candidates: {disagreement_total:,}
- Core common-hard high-experiment/low-all-core cases: {len(common_high):,}

## Decision and limitations

The gate remains `WAIT - BENCHMARK STILL INCOMPLETE`. The main blockers are
the empty assay-aligned construct benchmark and the absence of a runnable
Linux-compatible DeepPBS execution. FrozenPLM did not improve designed ranking
over SimplePC or the best sequence-only proxy.
"""
    (V042_RESULTS / "V0_4_2_PROGRESS.md").write_text(progress, encoding="utf-8")

    print(gate_text)
    print(report)


if __name__ == "__main__":
    main()
