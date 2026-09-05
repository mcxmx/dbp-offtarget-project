from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import project_root
from src.v0_5_phase4_diagnostics import (
    add_target_controls,
    assay_reference_control_sensitivity,
    build_hard_case_model_summary,
    build_hard_case_tables,
    build_wide_predictions,
    embedding_similarity_table,
    embedding_pca_table,
    load_reference_hard_cases,
    m0_shared_signal_correlations,
    partial_spearman_target_diagnostic,
    replay_primary_predictions,
    replay_validation,
    target_similarity_bin_performance,
    target_similarity_correlations,
    write_primary_frozen_manifest,
)
from src.v0_5_training import V05Config, load_v05_data


ROOT = project_root()
RESULTS = ROOT / "results" / "v0_5"
HARD_CASES = RESULTS / "hard_cases"
DOCS = ROOT / "docs" / "v0_5"


def git_revision(argument: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--verify", argument],
            cwd=ROOT,
            text=True,
        ).strip()
    except subprocess.CalledProcessError:
        return "unknown"


def fmt(value: object, digits: int = 4) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}f}"
    return str(value)


def markdown_table(frame: pd.DataFrame, columns: list[str] | None = None, digits: int = 4) -> str:
    table = frame if columns is None else frame[columns]
    lines = [
        "| " + " | ".join(str(column) for column in table.columns) + " |",
        "| " + " | ".join("---" for _ in table.columns) + " |",
    ]
    for _, row in table.iterrows():
        lines.append("| " + " | ".join(fmt(value, digits) for value in row) + " |")
    return "\n".join(lines)


def make_target_control_summary(wide: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein, group in wide.groupby("dbp_id", sort=True):
        control_values = {
            "TargetHamming": group["target_hamming"],
            "TargetEdit": group["target_edit"],
            "TargetKmerOverlap": group["target_kmer_overlap"],
        }
        correlations = {
            name: _safe_spearman(group["experimental_score"], values)
            for name, values in control_values.items()
        }
        rows.append(
            {
                "dbp_id": protein,
                "TargetHamming": correlations["TargetHamming"],
                "TargetEdit": correlations["TargetEdit"],
                "TargetKmerOverlap": correlations["TargetKmerOverlap"],
                "target_control_best": max(correlations.values()),
                "M2": _safe_spearman(group["experimental_score"], group["M2_score"]),
                "M3": _safe_spearman(group["experimental_score"], group["M3_score"]),
            }
        )
    return pd.DataFrame(rows)


def _safe_spearman(left: object, right: object) -> float:
    from src.v0_5_phase4_diagnostics import safe_spearman

    return safe_spearman(left, right)


def write_phase4_docs(
    *,
    hard_summary: pd.DataFrame,
    hard_by_protein: pd.DataFrame,
    subset_summary: pd.DataFrame,
    target_corr: pd.DataFrame,
    target_bins: pd.DataFrame,
    residual: pd.DataFrame,
    shuffle_p: pd.DataFrame,
    shuffle_t: pd.DataFrame,
    train_test: pd.DataFrame,
    pair_coverage: pd.DataFrame,
    embedding_similarity: pd.DataFrame,
    m0_shared: pd.DataFrame,
    target_controls: pd.DataFrame,
    assay_controls: pd.DataFrame,
    replay_check: pd.DataFrame,
) -> str:
    hard_total = int(hard_summary["reference_total"].iloc[0])
    subset_counts = subset_summary.set_index("subset")["count"].to_dict()
    joint_count = int(subset_counts.get("joint_controls_fail_m3_success", 0))
    represented_joint = subset_summary.loc[
        subset_summary["subset"].eq("joint_controls_fail_m3_success"), "proteins_represented"
    ].iloc[0] if joint_count else "none"
    all_fail_count = int(subset_counts.get("all_current_models_fail", 0))
    m3_hard = hard_summary.loc[hard_summary["model"].eq("M3")].iloc[0]
    m3_median_shuffle_corr = float(
        shuffle_p.loc[shuffle_p["condition"].eq("shuffled_protein"), "prediction_correlation_to_original"].median()
    )
    m3_median_target_shuffle_corr = float(
        shuffle_t.loc[shuffle_t["condition"].eq("shuffled_target"), "prediction_correlation_to_original"].median()
    )
    m3_target_corr = target_corr.loc[target_corr["model"].eq("M3"), "spearman_prediction_vs_target_kmer_overlap"]
    m2_target_corr = target_corr.loc[target_corr["model"].eq("M2"), "spearman_prediction_vs_target_kmer_overlap"]
    m1_target_corr = target_corr.loc[target_corr["model"].eq("M1"), "spearman_prediction_vs_target_kmer_overlap"]
    m3_residual = residual.loc[
        residual["model"].eq("M3"), "partial_spearman_residual_association"
    ]
    train_m3 = train_test.query("model == 'M3' and partition == 'train'")["spearman"]
    test_m3 = train_test.query("model == 'M3' and partition == 'test'")["spearman"]
    candidate_coverage = pair_coverage["candidate_coverage_fraction"]
    rank_coverage = pair_coverage["n_decile_bins_covered"]
    replay_max_error = float(replay_check["absolute_difference"].max())

    gate = "MODIFY"
    gate_reason = (
        "M3 does not beat the capacity-matched M1c or the target-blind M2 in the "
        "primary paired result, and the hard-case analysis does not provide sufficient "
        "evidence for a stable target-conditioned advantage. This rejects the current "
        "minimal implementation as the next training target, not the entire P,T,D concept."
    )
    if (
        joint_count == 0
        and m3_median_shuffle_corr >= 0.995
        and m3_median_target_shuffle_corr >= 0.995
    ):
        gate = "NO-GO"
        gate_reason = (
            "The frozen M3 produces no joint-control unique wins and its predictions are "
            "insensitive to both protein and target shuffles. The current target-conditioned "
            "concept has no empirical support under this benchmark."
        )

    hard_table = hard_summary[["model", "eligible", "resolved", "unresolved", "resolution_rate"]].copy()
    subset_table = subset_summary[
        [
            "subset",
            "count",
            "fraction_of_reference",
            "n_proteins",
            "proteins_represented",
            "experimental_E_score_median",
            "target_kmer_overlap_median",
        ]
    ].copy()
    target_medians = (
        target_corr.groupby("model", as_index=False)[
            ["spearman_prediction_vs_target_kmer_overlap", "spearman_prediction_vs_experimental"]
        ]
        .median()
        .rename(
            columns={
                "spearman_prediction_vs_target_kmer_overlap": "median_vs_target_kmer",
                "spearman_prediction_vs_experimental": "median_vs_experimental",
            }
        )
    )
    shuffle_summary = pd.DataFrame(
        [
            {
                "diagnostic": "shuffled protein",
                "median_prediction_correlation": m3_median_shuffle_corr,
                "median_experimental_spearman": float(
                    shuffle_p.loc[shuffle_p["condition"].eq("shuffled_protein"), "experimental_spearman"].median()
                ),
                "retrained": False,
            },
            {
                "diagnostic": "shuffled target",
                "median_prediction_correlation": m3_median_target_shuffle_corr,
                "median_experimental_spearman": float(
                    shuffle_t.loc[shuffle_t["condition"].eq("shuffled_target"), "experimental_spearman"].median()
                ),
                "retrained": False,
            },
        ]
    )
    train_test_summary = (
        train_test.groupby(["model", "partition"], as_index=False)["spearman"].median()
        .pivot(index="model", columns="partition", values="spearman")
        .reset_index()
        .rename(columns={"train": "train_median_spearman", "test": "test_median_spearman"})
    )
    embedding_diag = embedding_similarity.loc[
        ~embedding_similarity["same_dbp"],
        ["dbp_id_i", "dbp_id_j", "embedding_cosine_similarity", "embedding_euclidean_distance"],
    ].copy()
    protein_cluster_count = 4
    target_group_count = 4

    failure_doc = f"""# v0.5 Failure Analysis

## Scope and freeze

This is a secondary diagnostic over the frozen v0.5 primary evaluation. The
primary result files are protected by
`results/v0_5/PRIMARY_RESULTS_FROZEN_MANIFEST.txt` and were not overwritten.
Candidate-level scores were obtained by an exact replay of the frozen
fold/seed/config training protocol because Phase 3 did not persist predictions.
The replay was validated against the frozen seed-level Spearman values; the
maximum absolute difference was `{replay_max_error:.3e}`.

The hard-case reference is exactly the existing v0.3.1 set of **{hard_total:,}**
sequence-vs-experiment disagreement candidates. No new candidate threshold was
chosen after inspecting M0-M3. Resolved means seed-mean prediction percentile
at least 0.90 within protein, matching the existing v0.4.2 rule.

## Hard-case resolution

{markdown_table(hard_table)}

Per-protein counts are in `results/v0_5/hard_case_resolution_by_protein.csv`.
The denominator is consistent across the five complete models: `{hard_total:,}`
eligible candidates per model.

## Unique M3 wins and common failures

{markdown_table(subset_table)}

`joint_controls_fail_m3_success` contains **{joint_count}** candidates from
`{represented_joint}`. `all_current_models_fail.csv` contains **{all_fail_count}**
candidates where M0, M1, M1c, M2, and M3 all remain unresolved under the
predeclared triage rule. These are ranking discrepancies, not biological
binding-failure claims.

The individual row-level files are:

- `m1_fail_m3_success.csv`
- `m1c_fail_m3_success.csv`
- `m2_fail_m3_success.csv`
- `joint_controls_fail_m3_success.csv`
- `all_current_models_fail.csv`

## Target-similarity artifact audit

{markdown_table(target_medians)}

The per-protein and bin-level tables are
`target_similarity_model_correlations.csv` and
`target_similarity_bin_performance.csv`. Target-similarity tertiles were fixed
per protein before reading outcome values. M3's median prediction/target-k-mer
correlation was `{float(m3_target_corr.median()):.4f}`, compared with
M2 `{float(m2_target_corr.median()):.4f}` and M1 `{float(m1_target_corr.median()):.4f}`.

The residualized M3 diagnostic is in
`results/v0_5/partial_spearman_target_diagnostic.csv` and has median residual
association `{float(m3_residual.median()):.4f}`. This is a secondary rank
residual diagnostic, not a replacement for the primary Spearman metric.

## Protein and target use

{markdown_table(shuffle_summary)}

The shuffled-protein and shuffled-target diagnostics use already trained M3
models and set `retrained=False`; no shuffle result was used for training,
checkpoint selection, or tuning. Median prediction correlation after replacing
protein P was `{m3_median_shuffle_corr:.4f}`. After replacing target T it was
`{m3_median_target_shuffle_corr:.4f}`. A near-one correlation means the frozen
model's output is largely insensitive to that input; a large change without
better experimental ranking means the input affects predictions but does not
generalize in the intended direction.

## Why M0 can be stronger

{markdown_table(train_test_summary)}

M0 has no cross-protein protein representation to overfit. The protein-aware
heads are trained on only the designed train proteins in each LOCO fold, while
the held-out protein is a new biological unit. Their median training/test
diagnostic values are in `train_vs_test_performance.csv`; training values are
descriptive only and were not used for selection.

Pair sampling used `{float(pair_coverage["pairs_sampled"].median()):.0f}` pairs
per training protein. The median fraction of the 8,192 candidates appearing in
at least one sampled pair was `{float(candidate_coverage.median()):.1%}`, and
the median number of covered experimental rank deciles was
`{float(rank_coverage.median()):.1f}/10`. This is a plausible data-efficiency
bottleneck, but this phase does not increase pair counts.

M0's shared-signal diagnostics are in
`results/v0_5/m0_shared_signal_correlations.csv`, including correlations with
GC fraction and the prior k-mer3 proxy.

## SimplePC versus current M1

The prior v0.4.1 SimplePC result and current v0.5 M1 are not matched
experiments. SimplePC was trained on natural PBM data and evaluated on designed
proteins externally, while current M1 is trained only on the designed
protein-cluster LOCO training proteins. They also differ in protein
representation, DNA features, objective, pair sampling, normalization, and
evaluation regime. Therefore the apparent `{0.3616:.4f}` versus current
`{float(pd.read_csv(RESULTS / "primary_macro_summary.csv").query("model == 'M1'")["all7_macro_median"].iloc[0]):.4f}`
gap is a protocol comparison, not evidence that one implementation is
intrinsically better.

No bridge experiment was run: it would require a new training comparison and
could be mistaken for an optimization after the frozen primary result.

## Decision

**{gate}**

{gate_reason}

This gate is about the current frozen minimal model family and evidence needed
for the next experiment. It does not establish that all target-conditioned
architectures are impossible.
"""

    decision_doc = f"""# v0.5 Decision Memo

## Gate

**{gate}**

The decision is based on the frozen primary evaluation plus the predeclared
v0.4.2 hard-case resolution rule. The full candidate-level replay matched the
primary seed metrics within `{replay_max_error:.3e}`.

## Evidence

- M3 joint capacity-controlled/target-only unique wins: **{joint_count}**
  candidates, represented proteins: `{represented_joint}`.
- All five current models unresolved on: **{all_fail_count}** hard-case candidates.
- M3 hard-case resolution: `{int(m3_hard["resolved"]):,}/{int(m3_hard["eligible"]):,}`
  (`{float(m3_hard["resolution_rate"]):.1%}`).
- Median M3 prediction correlation after protein shuffle:
  `{m3_median_shuffle_corr:.4f}`.
- Median M3 prediction correlation after target shuffle:
  `{m3_median_target_shuffle_corr:.4f}`.
- Median residualized M3 association after controlling TargetKmerOverlap:
  `{float(m3_residual.median()):.4f}`.

## One primary next hypothesis

**The frozen global ESM protein representation and low-capacity conditioning
head are the primary bottleneck to test before changing the target-conditioned
concept.**

Rationale: replacing P with a deterministic different designed-protein
embedding changed the frozen M3 prediction by a median correlation of
`{m3_median_shuffle_corr:.4f}`, while replacing T changed it by
`{m3_median_target_shuffle_corr:.4f}`. The protein-aware heads also show a
training/test gap in `train_vs_test_performance.csv`, and M3 does not beat the
capacity-matched M1c or target-only M2 in the frozen primary result. Pair
coverage is a separate limitation: 512 pairs expose a median
`{float(candidate_coverage.median()):.1%}` of candidates and all ten rank
deciles, so it is recorded but is not selected as the primary next hypothesis.

Minimal falsification experiment: keep the target manifest, primary LOCO split,
RC semantics, training objective, and evaluation fixed; replace the single
global mean-pooled protein vector with a pre-registered residue-level or
protein-DNA local representation, without using designed test outcomes for
representation selection. Compare it against the frozen global-embedding M3
on the same folds and seeds.

Falsifier: if a protein representation with demonstrably non-constant
shuffled-P sensitivity still fails to improve M3 relative to M1c and M2 on
previously unseen proteins, and does not produce a stable capacity-controlled
hard-case advantage, then representation alone is not the primary bottleneck.

## Not concluded

The analysis does not prove that target conditioning is biologically invalid.
It shows that the current minimal FiLM model, frozen ESM global embedding,
ranking objective, and sparse pair protocol do not yet provide robust evidence
for a target-conditioned advantage.
"""

    simplepc_doc = """# SimplePC versus Current v0.5 M1 Audit

The prior v0.4.1 `SimpleProteinConditionalBaseline` and current v0.5 `M1`
numbers must not be interpreted as a head-to-head model comparison.

| Dimension | Prior SimplePC | Current v0.5 M1 |
| --- | --- | --- |
| Training proteins | Natural PBM benchmark | Designed DBP proteins in each LOCO training fold |
| Test regime | Designed external evaluation | Held-out designed protein cluster |
| Protein input | Frozen/compressed composition-style representation | Frozen ESM-2 t12 35M, 480-dimensional embedding |
| DNA input | Prior simple sequence feature representation | RC-symmetric one-hot 7-mer |
| Objective | Prior v0.4.1 regression/ridge-style protocol | Within-protein logistic pairwise ranking |
| Pair sampling | Prior protocol | 512 deterministic within-protein pairs |
| Evaluation unit | Prior designed external set | Protein first, 8,192 RC classes per protein |
| Target input | Not a target-conditioned model | M1 does not receive T |
| Leakage controls | v0.4.1 natural/designed protocol | Protein-cluster LOCO, no row-level split |

The most important difference is training regime: the prior result transfers
from a larger natural PBM collection, whereas current M1 learns its head from
only the other designed proteins in each fold. The comparison therefore mixes
training distribution, protein holdout, loss, representation, and sample
efficiency. It cannot identify a single causal reason for the numerical gap.

The current Phase 4 analysis deliberately does not run a bridge experiment.
Any bridge would be a new training diagnostic and would not alter the frozen
primary result.
"""

    (DOCS / "V0_5_FAILURE_ANALYSIS.md").write_text(failure_doc, encoding="utf-8")
    (DOCS / "V0_5_DECISION_MEMO.md").write_text(decision_doc, encoding="utf-8")
    (DOCS / "SIMPLEPC_VS_M1_AUDIT.md").write_text(simplepc_doc, encoding="utf-8")
    return gate


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    HARD_CASES.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    config = V05Config()
    manifest = RESULTS / "PRIMARY_RESULTS_FROZEN_MANIFEST.txt"
    if not manifest.exists():
        write_primary_frozen_manifest(
            ROOT,
            source_commit=git_revision("HEAD"),
            primary_result_commit=git_revision("v0.5-primary-evaluation-freeze^{}"),
            tag="v0.5-primary-evaluation-freeze",
            seeds=config.evaluation_seeds,
        )

    benchmark, targets, embedding_map, _ = load_v05_data()
    replay_path = RESULTS / "phase4_primary_predictions.parquet"
    if replay_path.exists():
        replay = {
            "predictions": pd.read_parquet(replay_path),
            "health": pd.read_csv(RESULTS / "phase4_replay_training_health.csv"),
            "train_test": pd.read_csv(RESULTS / "train_vs_test_performance.csv"),
            "pair_coverage": pd.read_csv(RESULTS / "pair_sampling_coverage.csv"),
            "shuffled_protein": pd.read_csv(RESULTS / "shuffled_protein_diagnostic.csv"),
            "shuffled_target": pd.read_csv(RESULTS / "shuffled_target_diagnostic.csv"),
        }
        if len(replay["predictions"]) != 8192 * 7 * 5 * len(config.evaluation_seeds):
            raise RuntimeError("Existing replay artifact has unexpected candidate coverage")
    else:
        replay = replay_primary_predictions(config)
        replay["predictions"].to_parquet(replay_path, index=False)
        replay["health"].to_csv(RESULTS / "phase4_replay_training_health.csv", index=False)
        replay["train_test"].to_csv(RESULTS / "train_vs_test_performance.csv", index=False)
        replay["pair_coverage"].to_csv(RESULTS / "pair_sampling_coverage.csv", index=False)
        replay["shuffled_protein"].to_csv(RESULTS / "shuffled_protein_diagnostic.csv", index=False)
        replay["shuffled_target"].to_csv(RESULTS / "shuffled_target_diagnostic.csv", index=False)
    replay_predictions = replay["predictions"]

    primary_seed = pd.read_csv(RESULTS / "primary_seed_level_results.csv")
    replay_check = replay_validation(
        replay_predictions,
        primary_seed,
        benchmark,
        RESULTS / "phase4_replay_validation.csv",
    )
    if not replay_check["matches_frozen_primary"].all():
        raise RuntimeError("Frozen primary replay did not reproduce one or more seed-level Spearman values")

    wide = add_target_controls(
        build_wide_predictions(replay_predictions, benchmark),
        targets,
    )
    wide.to_parquet(RESULTS / "phase4_wide_predictions.parquet", index=False)
    reference = load_reference_hard_cases(ROOT)
    hard = build_hard_case_tables(wide, reference, HARD_CASES)
    hard_summary = build_hard_case_model_summary(hard["hard"], RESULTS / "hard_case_model_summary.csv")
    target_corr = target_similarity_correlations(
        wide,
        RESULTS / "target_similarity_model_correlations.csv",
    )
    target_bins = target_similarity_bin_performance(
        wide,
        RESULTS / "target_similarity_bin_performance.csv",
    )
    residual = partial_spearman_target_diagnostic(
        wide,
        RESULTS / "partial_spearman_target_diagnostic.csv",
    )
    embedding_similarity = embedding_similarity_table(
        embedding_map,
        targets,
        RESULTS / "protein_embedding_similarity.csv",
    )
    embedding_pca_table(
        embedding_map,
        targets,
        RESULTS / "protein_embedding_pca.csv",
    )
    m0_shared = m0_shared_signal_correlations(
        wide,
        RESULTS / "m0_shared_signal_correlations.csv",
    )
    assay_controls = assay_reference_control_sensitivity(
        benchmark,
        targets,
        RESULTS / "assay_reference_control_sensitivity.csv",
    )
    target_controls = make_target_control_summary(wide)
    target_controls.to_csv(RESULTS / "target_control_summary.csv", index=False)
    gate = write_phase4_docs(
        hard_summary=hard_summary,
        hard_by_protein=hard["by_protein"],
        subset_summary=hard["summary"],
        target_corr=target_corr,
        target_bins=target_bins,
        residual=residual,
        shuffle_p=replay["shuffled_protein"],
        shuffle_t=replay["shuffled_target"],
        train_test=replay["train_test"],
        pair_coverage=replay["pair_coverage"],
        embedding_similarity=embedding_similarity,
        m0_shared=m0_shared,
        target_controls=target_controls,
        assay_controls=assay_controls,
        replay_check=replay_check,
    )
    print(f"Phase 4 diagnostics complete; gate={gate}")
    print(hard_summary.to_string(index=False))
    print(hard["summary"].to_string(index=False))
    print("replay max error:", replay_check["absolute_difference"].max())
    print("joint unique wins:", len(hard["joint_controls_fail_m3_success"]))


if __name__ == "__main__":
    main()
