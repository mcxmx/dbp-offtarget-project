from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import project_root


ROOT = project_root()
RESULTS = ROOT / "results" / "v0_5"
SEEDS = "17|29|43"


def fmt(value: object, digits: int = 4) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value))
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, int)):
        return f"{float(value):.{digits}f}"
    return str(value)


def markdown_table(frame: pd.DataFrame, columns: list[str] | None = None) -> str:
    table = frame if columns is None else frame[columns]
    headers = [str(column) for column in table.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in table.iterrows():
        values = []
        for value in row:
            values.append(fmt(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def metric_table(per_protein: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "dbp_id",
        "M0",
        "M1",
        "M1c",
        "M2",
        "M3",
        "delta_m3_minus_m1",
        "delta_m3_minus_m1c",
        "delta_m3_minus_m2",
    ]
    table = per_protein[columns].copy()
    table.columns = [
        "DBP",
        "M0",
        "M1",
        "M1c",
        "M2",
        "M3",
        "Delta M3-M1",
        "Delta M3-M1c",
        "Delta M3-M2",
    ]
    return table


def write_report() -> None:
    primary = pd.read_csv(RESULTS / "primary_per_protein_results.csv")
    strict = pd.read_csv(RESULTS / "strict_component_per_protein_results.csv")
    primary_macro = pd.read_csv(RESULTS / "primary_macro_summary.csv")
    strict_macro = pd.read_csv(RESULTS / "strict_component_macro_summary.csv")
    controls = pd.read_csv(RESULTS / "target_relative_controls_full.csv")
    context = pd.read_csv(RESULTS / "baseline_context_table.csv")
    health = pd.read_csv(RESULTS / "primary_training_health.csv")

    primary_all = metric_table(primary.sort_values("dbp_id"))
    primary_unseen = metric_table(
        primary.loc[primary["dbp_id"].isin(["DBP5", "DBP35", "DBP48", "DBP6", "DBP9"])]
        .sort_values("dbp_id")
    )
    strict_table = metric_table(strict.sort_values("dbp_id"))

    macro_columns = [
        "model",
        "all7_macro_median",
        "unseen5_macro_median",
        "proteins_evaluated",
        "unseen5_proteins_evaluated",
        "seeds",
    ]
    strict_macro_columns = [
        "model",
        "all_evaluated_macro_median",
        "all_evaluated_macro_mean",
        "all_evaluated_macro_sd",
        "proteins_evaluated",
        "components_evaluated",
        "seeds",
    ]

    control_pivot = controls.pivot(index="protein_id", columns="control", values="spearman").reset_index()
    control_pivot = control_pivot.rename(
        columns={
            "protein_id": "DBP",
            "TargetHamming": "Hamming",
            "TargetEdit": "Edit",
            "TargetKmerOverlap": "Kmer overlap",
        }
    ).sort_values("DBP")

    health_summary = pd.DataFrame(
        [
            {
                "quantity": "primary model-runs",
                "value": int((health["split_type"] == "protein_cluster_loco").sum()),
            },
            {
                "quantity": "strict model-runs",
                "value": int((health["split_type"] == "combined_component_loco").sum()),
            },
            {
                "quantity": "minimum prediction variance",
                "value": float(health["prediction_variance"].min()),
            },
            {
                "quantity": "maximum NaN/Inf count",
                "value": int(health["nan_inf_count"].max()),
            },
            {
                "quantity": "total training runtime (seconds)",
                "value": float(health["runtime_seconds"].sum()),
            },
            {
                "quantity": "training seeds",
                "value": SEEDS,
            },
        ]
    )

    m3 = primary_macro.loc[primary_macro["model"].eq("M3")].iloc[0]
    strict_m3 = strict_macro.loc[strict_macro["model"].eq("M3")].iloc[0]
    primary_m3_counts = {
        "M1": int(m3["m3_improved_over_m1_all7_count"]),
        "M1c": int(m3["m3_improved_over_m1c_all7_count"]),
        "M2": int(m3["m3_improved_over_m2_all7_count"]),
    }
    unseen_m3_counts = {
        "M1": int(m3["m3_improved_over_m1_unseen5_count"]),
        "M1c": int(m3["m3_improved_over_m1c_unseen5_count"]),
        "M2": int(m3["m3_improved_over_m2_unseen5_count"]),
    }

    lines = [
        "# v0.5 Primary Evaluation Results",
        "",
        "## Experimental status",
        "",
        "This report is the frozen v0.5 primary evidence artifact. It reports "
        "the complete four-fold protein-cluster LOCO evaluation and the three-fold "
        "combined-component sensitivity evaluation using the fixed seeds 17, 29, "
        "and 43.",
        "",
        "`protein_cluster_loco_fold_1` was previously used for engineering smoke "
        "training and is therefore marked `development_exposed` for DBP1/DBP3. "
        "Folds 2-4 cover the previously unseen five proteins DBP5, DBP35, DBP48, "
        "DBP6, and DBP9. No model architecture, target definition, split, "
        "hyperparameter, or pair protocol was changed after the smoke result.",
        "",
        "All metrics are calculated per protein on 8,192 canonical "
        "reverse-complement classes. Seed aggregation occurs within each "
        "protein/model before macro summarization. These results are not a "
        "row-level significance analysis.",
        "",
        "## Primary all-7",
        "",
        markdown_table(primary_all),
        "",
        "## Previously unseen five",
        "",
        markdown_table(primary_unseen),
        "",
        "## Primary macro summary",
        "",
        markdown_table(primary_macro, macro_columns),
        "",
        f"M3 improved over M1 on {primary_m3_counts['M1']}/7 proteins, over "
        f"M1c on {primary_m3_counts['M1c']}/7, and over M2 on "
        f"{primary_m3_counts['M2']}/7. For the previously unseen five, the "
        f"corresponding counts are {unseen_m3_counts['M1']}/5, "
        f"{unseen_m3_counts['M1c']}/5, and {unseen_m3_counts['M2']}/5.",
        "",
        "## Strict component sensitivity",
        "",
        "This is the assay-informed conservative sensitivity split. It controls "
        "the combined protein-cluster/target-group/motif leakage components and "
        "is not the primary deployment estimate.",
        "",
        markdown_table(strict_table),
        "",
        markdown_table(strict_macro, strict_macro_columns),
        "",
        f"Strict M3 median deltas were "
        f"{fmt(strict_m3['median_delta_m3_minus_m1'])} versus M1, "
        f"{fmt(strict_m3['median_delta_m3_minus_m1c'])} versus M1c, and "
        f"{fmt(strict_m3['median_delta_m3_minus_m2'])} versus M2. The direction "
        "is not uniformly preserved across proteins.",
        "",
        "## Target-relative controls",
        "",
        "These controls use only the independently sourced `primary_target`; "
        "they do not use PBM-derived motifs.",
        "",
        markdown_table(control_pivot),
        "",
        "## Seed stability",
        "",
        "The per-protein seed standard deviations are stored in "
        "`primary_per_protein_results.csv` and "
        "`strict_component_per_protein_results.csv`. The complete run health "
        "summary is:",
        "",
        markdown_table(health_summary),
        "",
        "## Baseline context",
        "",
        "Prior baselines are shown for context only. Coverage and training regime "
        "differ, so their macro values are not unconditional rankings.",
        "",
        markdown_table(context),
        "",
        "## Interpretation",
        "",
        "The frozen primary result does not provide robust evidence that M3 "
        "improves target-conditioned ranking over both protein-only matched "
        "controls. The unpaired all-7 macro medians are M1=0.0558, "
        "M1c=0.0605, and M3=0.0682, but the pre-registered paired "
        "per-protein median deltas are M3-M1=0.0124 and M3-M1c=-0.0227, "
        "with improvements on only 4/7 and 3/7 proteins respectively. The "
        "paired delta is the relevant comparison because the median of paired "
        "differences is not the difference of macro medians. M3-M2 has a "
        "paired median delta of -0.0201 and improves on only 3/7 proteins. "
        "The unseen-five subset shows small positive paired median deltas "
        "versus M1 and M1c but a negative delta versus M2, and is based on "
        "only five proteins.",
        "",
        "This report does not issue a final GO/NO-GO decision. Hard-case analysis "
        "and failure-resolution analysis are intentionally deferred to the next "
        "phase. The current result is a falsification-oriented benchmark of the "
        "frozen minimal model family, not evidence that a final proposed model "
        "should be implemented.",
        "",
        "## Future hypotheses",
        "",
        "- Analyze whether the limited M3 effect is concentrated in specific "
        "target/motif groups or sequence-distance regimes.",
        "- Separate implementation/capacity effects from genuine target-dependent "
        "ranking effects in the planned hard-case analysis.",
        "- Keep the strict component result as a leakage-sensitivity reference "
        "rather than replacing the primary LOCO result.",
        "",
    ]
    (RESULTS / "V0_5_PRIMARY_RESULTS.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    write_report()
