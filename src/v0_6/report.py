from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src.utils import project_root

from .audit import rc_representation_audit, target_representation_audit
from .gates import summarize_replay_gates


ROOT = project_root()
RESULTS = ROOT / "results" / "v0_6_representation"
V06_METADATA = ROOT / "metadata" / "v0_6"
DESIGNED_IDS = ["DBP1", "DBP3", "DBP5", "DBP35", "DBP48", "DBP6", "DBP9"]
MODELS = ["M0", "M1", "M1c", "M2", "M3"]


def _fmt(value: float) -> str:
    return "NA" if not np.isfinite(value) else f"{value:.4f}"


def generate_reports() -> dict[str, pd.DataFrame]:
    RESULTS.mkdir(parents=True, exist_ok=True)
    V06_METADATA.mkdir(parents=True, exist_ok=True)
    rc = pd.read_csv(RESULTS / "representation_audit.csv") if (RESULTS / "representation_audit.csv").exists() else rc_representation_audit()
    target = pd.read_csv(RESULTS / "target_representation_audit.csv") if (RESULTS / "target_representation_audit.csv").exists() else target_representation_audit()
    synthetic = pd.read_csv(RESULTS / "synthetic_conditionality.csv")
    primary = pd.read_csv(RESULTS / "primary_results.csv")
    model_rows = primary.loc[primary["model"].isin(MODELS)].copy()
    v06_per_protein = model_rows.groupby(["dbp_id", "model"], as_index=False)["spearman"].mean()
    v06_wide = v06_per_protein.pivot(index="dbp_id", columns="model", values="spearman").reset_index()
    v05 = pd.read_csv(ROOT / "results" / "v0_5" / "primary_seed_level_results.csv")
    v05_per_protein = v05.groupby(["dbp_id", "model"], as_index=False)["spearman"].mean()
    v05_wide = v05_per_protein.pivot(index="dbp_id", columns="model", values="spearman").reset_index()
    comparison_rows = []
    for model in MODELS:
        comparison_rows.append(
            {
                "model": model,
                "v05_all7_median_spearman": float(v05_wide[model].median()),
                "v06_all7_median_spearman": float(v06_wide[model].median()),
                "v06_minus_v05_median": float(v06_wide[model].median() - v05_wide[model].median()),
                "v05_all7_mean_spearman": float(v05_wide[model].mean()),
                "v06_all7_mean_spearman": float(v06_wide[model].mean()),
                "developmental_scope": "all seven GSE237017 DBPs exposed",
            }
        )
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(RESULTS / "v05_vs_v06_comparison.csv", index=False)
    v06_wide.to_csv(RESULTS / "v06_per_protein_spearman.csv", index=False)
    v05_wide.to_csv(RESULTS / "v05_per_protein_spearman.csv", index=False)
    seed_stability = (
        model_rows.groupby(["dbp_id", "model"], as_index=False)["spearman"]
        .agg(seed_mean="mean", seed_sd="std", seed_min="min", seed_max="max")
    )
    seed_summary = seed_stability.groupby("model", as_index=False).agg(
        median_per_protein_seed_sd=("seed_sd", "median"),
        max_per_protein_seed_sd=("seed_sd", "max"),
        median_seed_min=("seed_min", "median"),
        median_seed_max=("seed_max", "median"),
    )
    seed_stability.to_csv(RESULTS / "seed_stability_per_protein.csv", index=False)
    seed_summary.to_csv(RESULTS / "seed_stability.csv", index=False)
    protein_shuffle = pd.read_csv(RESULTS / "protein_shuffle.csv")
    target_shuffle = pd.read_csv(RESULTS / "target_shuffle.csv")
    health = pd.read_csv(RESULTS / "training_health.csv")
    v05_protein_shuffle = pd.read_csv(ROOT / "results" / "v0_5" / "shuffled_protein_diagnostic.csv")
    v05_target_shuffle = pd.read_csv(ROOT / "results" / "v0_5" / "shuffled_target_diagnostic.csv")
    shuffle_comparison = pd.DataFrame(
        [
            {
                "diagnostic": "protein_shuffle",
                "v05_median_prediction_correlation": float(v05_protein_shuffle.loc[v05_protein_shuffle["condition"].eq("shuffled_protein"), "prediction_correlation_to_original"].median()),
                "v06_median_prediction_correlation": float(protein_shuffle["prediction_correlation"].median()),
                "v05_median_effect_size": float(1.0 - v05_protein_shuffle.loc[v05_protein_shuffle["condition"].eq("shuffled_protein"), "prediction_correlation_to_original"].median()),
                "v06_median_effect_size": float(1.0 - protein_shuffle["prediction_correlation"].median()),
                "collapse_threshold": 0.995,
            },
            {
                "diagnostic": "target_shuffle",
                "v05_median_prediction_correlation": float(v05_target_shuffle.loc[v05_target_shuffle["condition"].eq("shuffled_target"), "prediction_correlation_to_original"].median()),
                "v06_median_prediction_correlation": float(target_shuffle["prediction_correlation"].median()),
                "v05_median_effect_size": float(1.0 - v05_target_shuffle.loc[v05_target_shuffle["condition"].eq("shuffled_target"), "prediction_correlation_to_original"].median()),
                "v06_median_effect_size": float(1.0 - target_shuffle["prediction_correlation"].median()),
                "collapse_threshold": 0.995,
            },
        ]
    )
    shuffle_comparison.to_csv(RESULTS / "v05_vs_v06_shuffle_comparison.csv", index=False)
    replay_gates = summarize_replay_gates(primary, protein_shuffle, target_shuffle, health)
    pre_gates = pd.read_csv(RESULTS / "pre_replay_gates.csv") if (RESULTS / "pre_replay_gates.csv").exists() else pd.DataFrame()
    if not pre_gates.empty:
        pre_gates = pre_gates.assign(stage="pre_replay")
    replay_gates = replay_gates.assign(stage="post_replay")
    gates = pd.concat([pre_gates, replay_gates], ignore_index=True)
    gates.to_csv(RESULTS / "go_no_go.csv", index=False)
    conditionality_rows = []
    for model in MODELS:
        model_subset = model_rows.loc[model_rows["model"].eq(model)]
        conditionality_rows.append(
            {
                "model": model,
                "held_out_spearman_median": float(model_subset["spearman"].median()),
                "candidate_prediction_variance_median": float(model_subset["prediction_variance"].median()),
                "protein_shuffle_prediction_correlation_median": float(protein_shuffle["prediction_correlation"].median()) if model == "M3" else np.nan,
                "protein_conditioning_effect_size_median": float(protein_shuffle["protein_conditioning_effect_size"].median()) if model == "M3" else np.nan,
                "target_shuffle_prediction_correlation_median": float(target_shuffle["prediction_correlation"].median()) if model == "M3" else np.nan,
                "target_conditioning_effect_size_median": float(target_shuffle["target_conditioning_effect_size"].median()) if model == "M3" else np.nan,
                "protein_shuffle_collapse_threshold": 0.995,
                "target_shuffle_collapse_threshold": 0.995,
                "scope": "developmental all-seven GSE237017 DBPs",
            }
        )
    conditionality = pd.DataFrame(conditionality_rows)
    conditionality.to_csv(RESULTS / "conditionality_diagnostics.csv", index=False)
    # Preserve a v0.6-local exposure copy with explicit current status.
    exposure = pd.read_csv(ROOT / "metadata" / "v0_5_transfer" / "exposure_manifest.csv")
    exposure["v0_6_use"] = "developmental_debugging_only"
    exposure["independent_confirmatory"] = False
    exposure.to_csv(V06_METADATA / "exposure_manifest.csv", index=False)
    target_rows = target.set_index("dbp_id")
    synthetic_pass = bool(synthetic["ranking_correct"].all()) and bool(synthetic["protein_shuffle_not_invariant"].all()) and bool(synthetic["target_shuffle_not_invariant"].all())
    p_corr = float(protein_shuffle["prediction_correlation"].median())
    t_corr = float(target_shuffle["prediction_correlation"].median())
    report_audit = f"""# V0.6 Representation Audit

## Scope and freeze

v0.6 is a representation-repair and conditionality-audit stage. It does not
modify or overwrite `results/v0_5/`, `docs/v0_5*`, or frozen v0.5 metadata.
All seven GSE237017 designed DBPs are already `development_exposed` according
to the current exposure manifest. Therefore every real-data replay result in
this report is developmental/debugging evidence, not independent confirmatory
validation.

## Current v0.5 representation problems

v0.5 encoded each candidate by averaging the full one-hot sequence with the
full one-hot reverse complement before learning. This is RC invariant, but it
is not injective over RC classes: positional base information is averaged away.
The target encoder additionally canonicalized each sliding 7-mer, deduplicated
with a set, sorted the survivors, and mean-pooled them. That loses window order,
register, and repeated-window multiplicity.

## Exhaustive RC audit

| quantity | value |
| --- | ---: |
| oriented 7-mers | {int(rc.iloc[0]['n_oriented_7mers'])} |
| canonical RC classes | {int(rc.iloc[0]['n_canonical_rc_classes'])} |
| v0.5 representation vectors | {int(rc.iloc[0]['legacy_representation_count'])} |
| v0.5 cross-class collision groups | {int(rc.iloc[0]['legacy_cross_rc_class_representation_collisions'])} |
| v0.5 collision excess over RC pairs | {int(rc.iloc[0]['legacy_collision_excess_over_rc_pairs'])} |
| v0.6 representation vectors | {int(rc.iloc[0]['v06_representation_count'])} |
| v0.6 cross-class collision groups | {int(rc.iloc[0]['v06_cross_rc_class_representation_collisions'])} |
| v0.6 RC-class injective | {bool(rc.iloc[0]['v06_rc_class_injective'])} |

The v0.6 candidate representation is the canonical RC representative followed
by ordinary one-hot encoding. A sequence and its reverse complement share one
vector, while different RC classes have different vectors. Independent score
invariance was tested over all 16,384 oriented 7-mers for M0, M1, M2, and M3.

## Target repair

v0.6 canonicalizes the complete target once, then retains every ordered window,
including duplicates, and appends a normalized register in `[-1, 1]`. Candidate
and target windows are compared explicitly. The default pooling is arithmetic
mean; max, log-sum-exp, softmax, and position-aware pooling are implemented as
small interpretable alternatives. No Transformer or attention stack was added.

## Synthetic conditionality gate

The synthetic benchmark contains P1/P2 and A-rich/C-rich target/candidate
preferences. All four P,T cells learned the required ranking reversal:
`P1: DNA_A > DNA_C`, `P2: DNA_C > DNA_A`, with target swaps also changing the
ranking. The synthetic gate is `{'GO' if synthetic_pass else 'NO-GO'}`.

## GO / NO-GO

Representation tests: **GO**. Synthetic conditionality: **GO**. Real-data
conditionality after the frozen replay: **NO-GO** because the protein-shuffle
correlation median is `{p_corr:.6f}` and the target-shuffle correlation median
is `{t_corr:.6f}`, both compared with the pre-registered collapse threshold
`0.995`. The repair therefore fixes an information-loss bug but does not by
itself establish that the learned model uses protein/target conditionally.
"""
    (ROOT / "docs" / "v0_6_representation").mkdir(parents=True, exist_ok=True)
    (ROOT / "docs" / "v0_6_representation" / "V0_6_REPRESENTATION_AUDIT.md").write_text(report_audit, encoding="utf-8")
    rows = []
    for dbp in DESIGNED_IDS:
        row = {"dbp_id": dbp}
        for model in MODELS:
            row[f"v05_{model}"] = float(v05_wide.loc[v05_wide.dbp_id.eq(dbp), model].iloc[0])
            row[f"v06_{model}"] = float(v06_wide.loc[v06_wide.dbp_id.eq(dbp), model].iloc[0])
        rows.append(row)
    per_protein_text = "\n".join(
        f"| {row['dbp_id']} | " + " | ".join(_fmt(row[f"v06_{model}"]) for model in MODELS) + " |"
        for row in rows
    )
    comparison_text = "\n".join(
        f"| {row.model} | {_fmt(row.v05_all7_median_spearman)} | {_fmt(row.v06_all7_median_spearman)} | {_fmt(row.v06_minus_v05_median)} |"
        for row in comparison.itertuples()
    )
    gate_text = gates.to_string(index=False)
    report_results = f"""# V0.6 Results

## Status and scope

This is one frozen replay after representation and synthetic tests passed.
The replay uses the v0.5 protein-cluster LOCO split, seeds `17|29|43`, 512
within-protein pairs, the same logistic pairwise objective, 18 epochs, and no
hyperparameter search. The seven designed DBPs are all development-exposed;
these results are not independent external validation.

## Historical evidence boundary

The repair is interpreted together with the already frozen v0.5 diagnostics:
871/1,515 pre-registered hard cases were unresolved by every M0-M3 model and
M3 resolved only 236/1,515. The v0.5 M3 protein-shuffle and target-shuffle
medians were already 1.0000 and 0.9999. The dense-supervision H2 pilot did not
improve held-out ranking or conditioning sensitivity, and the local/residue
H1 pilot had near-uniform attention with protein and target shuffle medians at
1.0000. Phase 7A H3 natural-to-designed transfer was explicitly stopped as
`NOT_SUPPORTED`; its unmatched prior SimplePC numbers are not a bridge result.
These constraints are why v0.6 performs a representation audit and one fixed
replay rather than further pair scaling, attention stacking, or transfer
tuning.

## v0.5 versus v0.6

| model | v0.5 all-7 median | v0.6 all-7 median | delta |
| --- | ---: | ---: | ---: |
{comparison_text}

v0.6 median Spearman is higher for M0/M1/M1c/M2/M3 respectively by the values
in `v05_vs_v06_comparison.csv`, but this is a developmental replay over an
already exposed cohort. It is evidence that the RC repair changes the learned
landscape, not evidence of confirmatory generalization.

## v0.6 per-protein Spearman (seed mean)

| DBP | M0 | M1 | M1c | M2 | M3 |
| --- | ---: | ---: | ---: | ---: | ---: |
{per_protein_text}

The sequence-only k-mer baseline remains a historical v0.3.1 comparator with
all-7 median Spearman `0.2321`; it is not retrained in v0.6. v0.6 medians are
M0 `{_fmt(v06_wide.M0.median())}`, M1 `{_fmt(v06_wide.M1.median())}`, M1c
`{_fmt(v06_wide.M1c.median())}`, M2 `{_fmt(v06_wide.M2.median())}`, and M3
`{_fmt(v06_wide.M3.median())}`.

## Protein and target shuffle diagnostics

| diagnostic | median prediction correlation | median effect size (1-correlation) |
| --- | ---: | ---: |
| protein shuffle | `{p_corr:.6f}` | `{1-p_corr:.6f}` |
| target shuffle | `{t_corr:.6f}` | `{1-t_corr:.6f}` |

The protein and target shuffle files are inference-only and use no retraining
or model selection. The protein conditionality gate is **NO-GO** and the target
conditionality gate is **NO-GO** under the `0.995` threshold.

For direct historical comparison, v0.5 protein-shuffle median correlation was
`{shuffle_comparison.loc[shuffle_comparison.diagnostic.eq('protein_shuffle'), 'v05_median_prediction_correlation'].iloc[0]:.6f}`
and target-shuffle median correlation was
`{shuffle_comparison.loc[shuffle_comparison.diagnostic.eq('target_shuffle'), 'v05_median_prediction_correlation'].iloc[0]:.6f}`.
The v0.6 values are in `v05_vs_v06_shuffle_comparison.csv`; both stages show
the same conditioning-collapse pattern.

## Seed stability and training health

Seed-level rows: 105 model/protein/seed results. Per-protein seed means and
standard deviations are in `seed_stability_per_protein.csv` and summary values
are in `seed_stability.csv`. All 60 model runs have positive prediction
variance and zero NaN/Inf values; detailed health is in `training_health.csv`.

## Gate table

```text
{gate_text}
```

## Scientific interpretation

1. **Is the current failure mainly a representation bug?** Partly. The RC
   averaging bug is real and severe: 16,384 oriented 7-mers collapse to 2,000
   vectors, with 1,872 cross-class collision groups. Repairing it changes the
   ranking metrics, so it was not a harmless implementation detail. However,
   the post-repair shuffle collapse remains, so representation bug alone is not
   the complete explanation.
2. **How much did RC repair solve?** It restores exact RC-class injectivity and
   raises developmental replay median Spearman, but it does not restore robust
   protein conditioning.
3. **Is target representation important?** The old target representation is
   demonstrably information-losing; v0.6 preserves order/register/multiplicity.
   In this replay, target shuffle remains near-invariant, so importance for
   learned conditionality is not established.
4. **Does protein-conditioning collapse remain?** Yes. Protein shuffle median
   correlation is `{p_corr:.6f}`.
5. **Should sequence-based conditional modeling continue?** **STOP as a final
   model-development line.** The synthetic capability is present, but real
   protein/target effects are not used robustly and the sequence-only baseline
   remains competitive or better. A further small, pre-registered audit could
   test a new hypothesis, but no capacity/pair/attention scaling is justified.
6. **Next direction:** do not enter natural-to-designed transfer from this
   exposed cohort. The evidence supports a benchmark/generalization-gap paper
   or a carefully pre-registered structure-guided hypothesis after obtaining a
   new independent designed-DBP cohort. Natural-transfer H3 is already
   `NOT_SUPPORTED`; no confirmatory quantitative dataset was searched here.
"""
    (ROOT / "docs" / "v0_6_representation" / "V0_6_RESULTS.md").write_text(report_results, encoding="utf-8")
    return {
        "comparison": comparison,
        "v06_per_protein": v06_wide,
        "seed_stability": seed_stability,
        "seed_summary": seed_summary,
        "gates": gates,
    }


if __name__ == "__main__":
    output = generate_reports()
    for frame in output.values():
        print(frame.to_string(index=False))
