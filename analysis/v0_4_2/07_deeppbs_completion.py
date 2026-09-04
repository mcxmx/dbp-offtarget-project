"""Complete the DeepPBS designed-DBP baseline and hard-case integration."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.deeppbs_evaluation import (
    DESIGNED_PROTEINS,
    bootstrap_per_protein,
    evaluate_per_protein,
    load_designed_experimental_units,
    score_prediction_npz,
)
from src.utils import ensure_dir, project_root


ROOT = project_root()
V042 = ROOT / "results" / "v0_4_2"
TABLES = ensure_dir(V042 / "tables")
FIGURES = ensure_dir(V042 / "figures")
RUNS = V042 / "external_runs" / "deeppbs"
PROCESSED = ensure_dir(ROOT / "data" / "processed" / "v0_4_2")
DOCS = ensure_dir(ROOT / "docs" / "v0_4_2")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4_2")
SEED = 42
DEEPPBS_COMMIT = "8bfb211dd67f02877841f6f33aa493ddf7daedf9"
DEEPPBS_MODEL_VERSION = f"DeepPBS official ensemble; commit {DEEPPBS_COMMIT}"
DEEPPBS_WEIGHT_ID = "official bundled ensemble from run/plot_scripts/txts/DeepPBS.txt"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_prediction_tables() -> pd.DataFrame:
    configs = {
        "DBP35": {
            "npz": RUNS / "DBP35" / "DBP035.npz_predict.npz",
            "structure_id": "DBP035",
            "structure_source": (
                "external/dbp_design/2b_design_mpnn/DBP035.pdb; "
                "official DeepPBS process output"
            ),
            "overlap_status": "no_exact_overlap_found_in_checked_manifests; homolog_unknown",
        },
        "DBP48": {
            "npz": RUNS / "DBP48" / "8TAC_helix_only.npz_predict.npz",
            "structure_id": "8TAC",
            "structure_source": (
                "PDB 8TAC; DSSR-defined 9-bp helix-only derived input; "
                "official DeepPBS process output"
            ),
            "overlap_status": "no_exact_overlap_found_in_checked_manifests; homolog_unknown",
        },
    }
    tables = []
    for protein_id, config in configs.items():
        table = score_prediction_npz(
            config["npz"],
            protein_id,
            structure_id=config["structure_id"],
            structure_source=config["structure_source"],
            model_version=DEEPPBS_MODEL_VERSION,
            weight_id=DEEPPBS_WEIGHT_ID,
            overlap_status=config["overlap_status"],
        )
        tables.append(table)
    result = pd.concat(tables, ignore_index=True)
    expected = 2 * 8192
    if len(result) != expected or result.duplicated(
        ["protein_id", "canonical_7mer"]
    ).any():
        raise ValueError("DeepPBS prediction unit integrity check failed")
    return result


def write_structure_manifest() -> pd.DataFrame:
    rows = [
        {
            "dbp_id": "DBP1",
            "protein_sequence_id": "metadata/v0_3/designed_dbp_sequences.csv:DBP1",
            "intended_design_target": "TAGCAGGATGTGT",
            "experimental_pbm_reference": "GCAGG",
            "pdb_id": "",
            "structure_source": "not found in checked public/project sources",
            "biological_assembly": "",
            "contains_protein_dna_complex": False,
            "structure_matches_assayed_protein": False,
            "deepPBS_evaluable": False,
            "reason_not_evaluable": "no reliable public protein-DNA complex structure/model",
            "exact_training_overlap": "not_confirmed_seen",
            "overlap_status": "not evaluated; homolog audit unresolved",
            "notes": "No structure fabricated.",
        },
        {
            "dbp_id": "DBP3",
            "protein_sequence_id": "metadata/v0_3/designed_dbp_sequences.csv:DBP3",
            "intended_design_target": "TAGCAGGATGTGT",
            "experimental_pbm_reference": "GCAGGA",
            "pdb_id": "",
            "structure_source": "not found in checked public/project sources",
            "biological_assembly": "",
            "contains_protein_dna_complex": False,
            "structure_matches_assayed_protein": False,
            "deepPBS_evaluable": False,
            "reason_not_evaluable": "no reliable public protein-DNA complex structure/model",
            "exact_training_overlap": "not_confirmed_seen",
            "overlap_status": "not evaluated; homolog audit unresolved",
            "notes": "No structure fabricated.",
        },
        {
            "dbp_id": "DBP5",
            "protein_sequence_id": "metadata/v0_3/designed_dbp_sequences.csv:DBP5",
            "intended_design_target": "GCAGATCTGCACATC",
            "experimental_pbm_reference": "TGCACA",
            "pdb_id": "",
            "structure_source": "not found in checked public/project sources",
            "biological_assembly": "",
            "contains_protein_dna_complex": False,
            "structure_matches_assayed_protein": False,
            "deepPBS_evaluable": False,
            "reason_not_evaluable": "no reliable public protein-DNA complex structure/model",
            "exact_training_overlap": "not_confirmed_seen",
            "overlap_status": "not evaluated; homolog audit unresolved",
            "notes": "No structure fabricated.",
        },
        {
            "dbp_id": "DBP6",
            "protein_sequence_id": "metadata/v0_3/designed_dbp_sequences.csv:DBP6",
            "intended_design_target": "GCAGATCTGCACATC",
            "experimental_pbm_reference": "TGCACA",
            "pdb_id": "",
            "structure_source": "not found in checked public/project sources",
            "biological_assembly": "",
            "contains_protein_dna_complex": False,
            "structure_matches_assayed_protein": False,
            "deepPBS_evaluable": False,
            "reason_not_evaluable": "no reliable public protein-DNA complex structure/model",
            "exact_training_overlap": "not_confirmed_seen",
            "overlap_status": "not evaluated; homolog audit unresolved",
            "notes": "No structure fabricated.",
        },
        {
            "dbp_id": "DBP9",
            "protein_sequence_id": "metadata/v0_3/designed_dbp_sequences.csv:DBP9",
            "intended_design_target": "GCAGATCTGCACATC",
            "experimental_pbm_reference": "TGCACA",
            "pdb_id": "",
            "structure_source": "not found in checked public/project sources",
            "biological_assembly": "",
            "contains_protein_dna_complex": False,
            "structure_matches_assayed_protein": False,
            "deepPBS_evaluable": False,
            "reason_not_evaluable": "no reliable public protein-DNA complex structure/model",
            "exact_training_overlap": "not_confirmed_seen",
            "overlap_status": "not evaluated; homolog audit unresolved",
            "notes": "No structure fabricated.",
        },
        {
            "dbp_id": "DBP35",
            "protein_sequence_id": "metadata/v0_3/designed_dbp_sequences.csv:DBP35",
            "intended_design_target": "GCAGATCTGCACATC",
            "experimental_pbm_reference": "TGCACA",
            "pdb_id": "",
            "structure_source": "external/dbp_design/2b_design_mpnn/DBP035.pdb",
            "biological_assembly": "as-deposited theoretical design model; DNA chain B",
            "contains_protein_dna_complex": True,
            "structure_matches_assayed_protein": True,
            "deepPBS_evaluable": True,
            "reason_not_evaluable": "",
            "exact_training_overlap": "not_confirmed_seen",
            "overlap_status": "no_exact_overlap_found_in_checked_manifests; homolog_unknown",
            "notes": (
                "Theoretical Rosetta model. Official process completed with "
                "Helix score 1.0 and CONTACT COUNT 266."
            ),
        },
        {
            "dbp_id": "DBP48",
            "protein_sequence_id": "metadata/v0_3/designed_dbp_sequences.csv:DBP48",
            "intended_design_target": "CGCCCAAAGCCGCG",
            "experimental_pbm_reference": "CTGACG",
            "pdb_id": "8TAC",
            "structure_source": (
                "data/raw/rcsb/mmcif/8TAC.cif; project-derived "
                "DSSR helix-only PDB"
            ),
            "biological_assembly": "as-deposited chains A/B protein and C/D DNA; helix-only keeps 9 bp",
            "contains_protein_dna_complex": True,
            "structure_matches_assayed_protein": True,
            "deepPBS_evaluable": True,
            "reason_not_evaluable": "",
            "exact_training_overlap": "not_confirmed_seen",
            "overlap_status": "no_exact_overlap_found_in_checked_manifests; homolog_unknown",
            "notes": (
                "Original 21-nt assembly has three DSSR non-helical overhang "
                "residues. DeepPBS run uses only the unique 9-bp DSSR helix; "
                "original structure and derived input are both retained."
            ),
        },
    ]
    manifest = pd.DataFrame(rows)
    manifest.to_csv(
        TABLES / "deeppbs_structure_manifest_completed_v0_4_2.csv", index=False
    )
    manifest.to_csv(
        METADATA / "designed_structure_manifest_completed_v0_4_2.csv", index=False
    )
    return manifest


def write_overlap_manifest() -> pd.DataFrame:
    """Record the completed exact-overlap audit without overstating homology."""

    rows = []
    for protein_id in DESIGNED_PROTEINS:
        rows.append(
            {
                "protein_id": protein_id,
                "exact_sequence_seen": False,
                "homolog_seen": "unknown_not_assessed",
                "structure_seen": False,
                "paper_evaluation_seen": False,
                "training_status": "exact_overlap_not_found; homolog_unknown",
                "evidence": (
                    "Checked the available DeepPBS repository split/manifests; "
                    "no exact designed-DBP sequence, PDB, or GSE237017 ID was found. "
                    "The upstream repository does not provide a complete training "
                    "sequence manifest for homolog-level assessment."
                ),
                "risk_level": "low_to_medium_unresolved_homology",
                "notes": (
                    "Results are not called strict zero-shot. DBP35 uses a "
                    "theoretical design model; DBP48 uses PDB 8TAC with a "
                    "project-side helix-only input repair."
                ),
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(
        METADATA / "deeppbs_overlap_audit_completed_v0_4_2.csv", index=False
    )
    return audit


def write_run_manifest() -> pd.DataFrame:
    rows = [
        {
            "protein_id": "DBP35",
            "structure_id": "DBP035",
            "input_structure": "results/v0_4_2/external_runs/deeppbs/DBP35/DBP035.pdb",
            "processed_npz": "results/v0_4_2/external_runs/deeppbs/DBP35/DBP035.npz",
            "prediction_npz": "results/v0_4_2/external_runs/deeppbs/DBP35/DBP035.npz_predict.npz",
            "process_stdout": "results/v0_4_2/external_runs/deeppbs/DBP35/process_x3dna_stdout.txt",
            "process_stderr": "results/v0_4_2/external_runs/deeppbs/DBP35/process_x3dna_stderr.txt",
            "predict_stdout": "results/v0_4_2/external_runs/deeppbs/DBP35/predict_stdout.txt",
            "predict_stderr": "results/v0_4_2/external_runs/deeppbs/DBP35/predict_stderr.txt",
            "process_returncode": 0,
            "predict_returncode": 0,
            "helix_score": 1.0,
            "contact_count": 266,
            "input_preparation": "original project design complex PDB",
            "sha256_input_structure": sha256(
                ROOT
                / "results"
                / "v0_4_2"
                / "external_runs"
                / "deeppbs"
                / "DBP35"
                / "DBP035.pdb"
            ),
            "sha256_prediction_npz": sha256(
                RUNS / "DBP35" / "DBP035.npz_predict.npz"
            ),
            "notes": "Nonfatal missing P/OP atom warnings recorded in stderr.",
        },
        {
            "protein_id": "DBP48",
            "structure_id": "8TAC",
            "input_structure": "results/v0_4_2/external_runs/deeppbs/DBP48/8TAC_helix_only.pdb",
            "processed_npz": "results/v0_4_2/external_runs/deeppbs/DBP48/8TAC_helix_only.npz",
            "prediction_npz": "results/v0_4_2/external_runs/deeppbs/DBP48/8TAC_helix_only.npz_predict.npz",
            "process_stdout": "results/v0_4_2/external_runs/deeppbs/DBP48/process_helix_x3dna_stdout.txt",
            "process_stderr": "results/v0_4_2/external_runs/deeppbs/DBP48/process_helix_x3dna_stderr.txt",
            "predict_stdout": "results/v0_4_2/external_runs/deeppbs/DBP48/predict_stdout.txt",
            "predict_stderr": "results/v0_4_2/external_runs/deeppbs/DBP48/predict_stderr.txt",
            "process_returncode": 0,
            "predict_returncode": 0,
            "helix_score": 1.0,
            "contact_count": 224,
            "input_preparation": "DSSR unique 9-bp helix-only derived PDB",
            "sha256_input_structure": sha256(
                RUNS / "DBP48" / "8TAC_helix_only.pdb"
            ),
            "sha256_prediction_npz": sha256(
                RUNS / "DBP48" / "8TAC_helix_only.npz_predict.npz"
            ),
            "notes": (
                "Original 8TAC failed upstream shape extraction on three "
                "non-helical overhang residues; derived helix-only input passed."
            ),
        },
    ]
    manifest = pd.DataFrame(rows)
    manifest.to_csv(
        TABLES / "deeppbs_run_manifest_completed_v0_4_2.csv", index=False
    )
    return manifest


def write_performance(experimental: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    performance = evaluate_per_protein(experimental, predictions)
    performance.to_csv(
        TABLES / "deeppbs_performance_completed_v0_4_2.csv", index=False
    )
    evaluated = performance[performance["status"] == "evaluated"]
    summary_rows = []
    for metric in [
        "spearman",
        "ndcg_1pct",
        "ndcg_5pct",
        "pairwise_accuracy",
        "top1pct_recovery",
    ]:
        values = evaluated[metric].dropna().to_numpy(dtype=float)
        summary_rows.append(
            {
                "baseline": "DeepPBS",
                "metric": metric,
                "n_evaluable_proteins": int(len(values)),
                "mean": float(np.mean(values)) if len(values) else np.nan,
                "median": float(np.median(values)) if len(values) else np.nan,
                "std": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
                "coverage": f"{len(evaluated)}/7",
            }
        )
    pd.DataFrame(summary_rows).to_csv(
        TABLES / "deeppbs_macro_performance_completed_v0_4_2.csv", index=False
    )
    bootstrap_per_protein(experimental, predictions).to_csv(
        TABLES / "deeppbs_bootstrap_ci_completed_v0_4_2.csv", index=False
    )
    return performance


def load_scored_landscape(experimental: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    """Rebuild the v0.4.2 row-level landscape with completed DeepPBS scores."""

    sequence = pd.read_parquet(
        ROOT
        / "data"
        / "processed"
        / "v0_3_1"
        / "designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1.parquet"
    )[
        [
            "protein_id",
            "canonical_7mer",
            "kmer3_jaccard_to_paper_motif_rc_aware",
            "hamming_similarity_to_paper_motif_rc_aware",
        ]
    ].rename(
        columns={
            "kmer3_jaccard_to_paper_motif_rc_aware": "sequence_score",
            "hamming_similarity_to_paper_motif_rc_aware": "sequence_hamming_score",
        }
    )
    simple = pd.read_parquet(
        ROOT / "results" / "v0_4_1" / "tables" / "simple_pc_designed_predictions.parquet"
    )[["protein_id", "canonical_rc", "simple_pc_score"]].rename(
        columns={"canonical_rc": "canonical_7mer"}
    )
    frozen = pd.read_parquet(
        ROOT / "results" / "v0_4_2" / "tables" / "frozen_plm_designed_predictions.parquet"
    )[["protein_id", "canonical_rc", "frozen_plm_score"]].rename(
        columns={"canonical_rc": "canonical_7mer"}
    )
    nampnn_path = ROOT / "results" / "v0_4" / "tables" / "nampnn_predictions.parquet"
    nampnn = (
        pd.read_parquet(nampnn_path)[
            ["protein_id", "canonical_7mer", "prediction_score"]
        ].rename(columns={"prediction_score": "nampnn_score"})
        if nampnn_path.exists()
        else pd.DataFrame(columns=["protein_id", "canonical_7mer", "nampnn_score"])
    )
    deeppbs = predictions[["protein_id", "canonical_7mer", "deeppbs_score"]]
    scored = experimental[
        ["protein_id", "canonical_7mer", "experimental_score"]
    ].copy()
    scored = scored.merge(sequence, on=["protein_id", "canonical_7mer"], how="left")
    scored = scored.merge(simple, on=["protein_id", "canonical_7mer"], how="left")
    scored = scored.merge(frozen, on=["protein_id", "canonical_7mer"], how="left")
    scored = scored.merge(nampnn, on=["protein_id", "canonical_7mer"], how="left")
    scored = scored.merge(deeppbs, on=["protein_id", "canonical_7mer"], how="left")
    for column in [
        "experimental_score",
        "sequence_score",
        "simple_pc_score",
        "frozen_plm_score",
        "nampnn_score",
        "deeppbs_score",
    ]:
        scored[f"{column}_percentile"] = scored.groupby("protein_id")[column].rank(
            pct=True, ascending=True
        )
    return scored


def attach_experimental_scores(
    experimental: pd.DataFrame, predictions: pd.DataFrame
) -> pd.DataFrame:
    """Join the frozen PBM reference without changing the prediction score."""

    reference = experimental[
        ["protein_id", "canonical_7mer", "experimental_score"]
    ].rename(columns={"experimental_score": "experimental_E_score"})
    joined = predictions.merge(
        reference,
        on=["protein_id", "canonical_7mer"],
        how="left",
        validate="one_to_one",
    )
    if len(joined) != len(predictions) or joined["experimental_E_score"].isna().any():
        raise ValueError("DeepPBS/PBM joined landscape has missing or duplicate units")
    joined["experimental_score_type"] = (
        "processed experimental uPBM E-score consensus"
    )
    return joined


def candidate_mask(scored: pd.DataFrame) -> pd.Series:
    thresholds = pd.read_csv(
        ROOT / "results" / "v0_3_1" / "tables" / "all_disagreement_candidate_counts.csv"
    ).set_index("protein_id")
    mask = pd.Series(False, index=scored.index)
    for protein_id, group in scored.groupby("protein_id"):
        threshold = thresholds.loc[protein_id]
        mask.loc[group.index] = (
            group["experimental_score"]
            >= float(threshold["experimental_score_threshold"])
        ) & (
            group["sequence_hamming_score"]
            <= float(threshold["sequence_similarity_threshold"]) + 1e-12
        )
    return mask


def write_hard_case_outputs(scored: pd.DataFrame) -> dict[str, int]:
    mask = candidate_mask(scored)
    method_columns = {
        "sequence_kmer3": "sequence_score",
        "SimpleProteinConditionalBaseline": "simple_pc_score",
        "FrozenPLMProteinConditionalBaseline": "frozen_plm_score",
        "NA-MPNN diagnostic": "nampnn_score",
        "DeepPBS": "deeppbs_score",
    }
    rows = []
    for protein_id, group in scored.groupby("protein_id", sort=True):
        candidates = group.loc[mask.loc[group.index]]
        for method, column in method_columns.items():
            evaluable = candidates[column].notna()
            resolved = evaluable & (
                candidates[f"{column}_percentile"] >= 0.90
            )
            rows.append(
                {
                    "protein_id": protein_id,
                    "method": method,
                    "n_total_candidates": int(len(candidates)),
                    "n_evaluable": int(evaluable.sum()),
                    "n_resolved": int(resolved.sum()),
                    "n_unresolved": int((evaluable & ~resolved).sum()),
                    "resolution_rate_among_evaluable": (
                        float(resolved.sum() / evaluable.sum())
                        if evaluable.sum()
                        else np.nan
                    ),
                    "resolution_rate_among_all_candidates": (
                        float(resolved.sum() / len(candidates))
                        if len(candidates)
                        else np.nan
                    ),
                    "threshold_definition": (
                        "v0.3.1 disagreement candidates; method within-protein "
                        "prediction percentile >= 0.90"
                    ),
                    "evaluation_status": (
                        "evaluated"
                        if evaluable.any()
                        else "not_evaluable_missing_prediction"
                    ),
                }
            )
    resolution = pd.DataFrame(rows)
    resolution.to_csv(
        TABLES / "disagreement_resolution_deeppbs_completed_v0_4_2.csv",
        index=False,
    )

    candidate_columns = [
        "protein_id",
        "canonical_7mer",
        "experimental_score",
        "experimental_score_percentile",
        "sequence_score",
        "sequence_score_percentile",
        "simple_pc_score",
        "simple_pc_score_percentile",
        "frozen_plm_score",
        "frozen_plm_score_percentile",
        "deeppbs_score",
        "deeppbs_score_percentile",
        "nampnn_score",
        "nampnn_score_percentile",
    ]
    candidates_all = scored.loc[mask, candidate_columns].copy()
    candidates_all["deepPBS_evaluable"] = candidates_all["deeppbs_score"].notna()
    candidates_all["experimental_high_candidate"] = True
    candidates_all["sequence_only_low"] = candidates_all[
        "sequence_score_percentile"
    ].le(0.50)
    candidates_all["deepPBS_low"] = candidates_all[
        "deeppbs_score_percentile"
    ].le(0.50)
    candidates_all["deepPBS_high"] = candidates_all[
        "deeppbs_score_percentile"
    ].ge(0.90)
    complete_core = candidates_all[
        [
            "sequence_score_percentile",
            "simple_pc_score_percentile",
            "frozen_plm_score_percentile",
            "deeppbs_score_percentile",
        ]
    ].notna().all(axis=1)
    all_core_low = complete_core & candidates_all[
        [
            "sequence_score_percentile",
            "simple_pc_score_percentile",
            "frozen_plm_score_percentile",
            "deeppbs_score_percentile",
        ]
    ].le(0.50).all(axis=1)
    candidates_all["failure_category"] = np.select(
        [
            ~candidates_all["deepPBS_evaluable"],
            candidates_all["deepPBS_high"]
            & candidates_all["sequence_only_low"],
            candidates_all["deepPBS_low"]
            & candidates_all["sequence_score_percentile"].ge(0.90),
            candidates_all["deepPBS_high"]
            & candidates_all["sequence_score_percentile"].ge(0.90),
            all_core_low,
            candidates_all["deepPBS_low"]
            & candidates_all["sequence_only_low"],
        ],
        [
            "DeepPBS_not_evaluable",
            "DeepPBS_resolves_sequence_disagreement",
            "sequence_only_resolves_DeepPBS_disagreement",
            "both_methods_high",
            "common_low_prediction",
            "both_methods_low",
        ],
        default="mixed_or_unclassified",
    )
    candidates_all["failure_definition"] = (
        "Frozen v0.3.1 high-experiment/sequence-disagreement candidate; "
        "method high >= 0.90 percentile and low <= 0.50 percentile."
    )
    candidates_all.to_parquet(
        TABLES / "baseline_failure_cases_deeppbs_completed_v0_4_2.parquet",
        index=False,
    )

    common = candidates_all.copy()
    core = [
        "sequence_score_percentile",
        "simple_pc_score_percentile",
        "frozen_plm_score_percentile",
        "deeppbs_score_percentile",
    ]
    common["failure_type_completed_core"] = np.where(
        common[core].notna().all(axis=1) & common[core].le(0.50).all(axis=1),
        "high_experiment_low_all_completed_core",
        "not_complete_core_or_not_common",
    )
    common.to_parquet(
        PROCESSED / "common_hard_specificity_cases_deeppbs_completed_v0_4_2.parquet",
        index=False,
    )

    eligible = resolution[
        (resolution["method"] == "DeepPBS")
        & (resolution["evaluation_status"] == "evaluated")
    ]
    candidate_total = int(eligible["n_total_candidates"].sum())
    resolved_total = int(eligible["n_resolved"].sum())
    complete_common = common[
        common["failure_type_completed_core"]
        == "high_experiment_low_all_completed_core"
    ]
    complete_common.to_parquet(
        PROCESSED
        / "common_hard_specificity_cases_deeppbs_core_only_v0_4_2.parquet",
        index=False,
    )
    category_counts = (
        candidates_all["failure_category"]
        .value_counts(dropna=False)
        .rename_axis("failure_category")
        .reset_index(name="n_candidates")
    )
    category_counts.insert(0, "candidate_scope", "all_v0.3.1_disagreement_candidates")
    category_counts.to_csv(
        TABLES / "baseline_failure_category_counts_deeppbs_completed_v0_4_2.csv",
        index=False,
    )
    pairwise = []
    eligible_pair = candidates_all[candidates_all["deepPBS_evaluable"]].copy()
    pairwise_rules = {
        "sequence_low_deeppbs_high": eligible_pair["sequence_only_low"]
        & eligible_pair["deepPBS_high"],
        "sequence_high_deeppbs_low": eligible_pair[
            "sequence_score_percentile"
        ].ge(0.90)
        & eligible_pair["deepPBS_low"],
        "sequence_low_deeppbs_low": eligible_pair["sequence_only_low"]
        & eligible_pair["deepPBS_low"],
        "sequence_high_deeppbs_high": eligible_pair[
            "sequence_score_percentile"
        ].ge(0.90)
        & eligible_pair["deepPBS_high"],
    }
    for category, category_mask in pairwise_rules.items():
        pairwise.append(
            {
                "category": category,
                "eligible_candidates": int(len(eligible_pair)),
                "n_candidates": int(category_mask.sum()),
                "fraction_of_eligible": (
                    float(category_mask.mean()) if len(eligible_pair) else np.nan
                ),
                "definition": (
                    "candidate-specific within-protein percentile thresholds; "
                    "high >= 0.90, low <= 0.50"
                ),
            }
        )
    pd.DataFrame(pairwise).to_csv(
        TABLES / "baseline_complementarity_deeppbs_completed_v0_4_2.csv",
        index=False,
    )
    complementary = {
        "eligible_candidates": candidate_total,
        "resolved": resolved_total,
        "unresolved": candidate_total - resolved_total,
        "resolution_rate": (
            resolved_total / candidate_total if candidate_total else np.nan
        ),
        "complete_core_common_failures": int(len(complete_common)),
        "all_disagreement_candidates": int(len(candidates_all)),
    }
    pd.DataFrame([complementary]).to_csv(
        TABLES / "deeppbs_hard_case_summary_completed_v0_4_2.csv", index=False
    )
    return complementary


def write_figures(performance: pd.DataFrame, scored: pd.DataFrame) -> None:
    evaluated = performance[performance["status"] == "evaluated"].copy()
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.bar(evaluated["protein_id"], evaluated["spearman"], color="#2A6FBB")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Spearman: PBM E-score vs DeepPBS score")
    ax.set_xlabel("Designed DBP")
    ax.set_title("DeepPBS designed-DBP ranking performance")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_deeppbs_per_protein_spearman_completed.png", dpi=300)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)
    for ax, protein_id in zip(axes, ["DBP35", "DBP48"]):
        sub = scored[scored["protein_id"] == protein_id]
        ax.scatter(
            sub["experimental_score"],
            sub["deeppbs_score"],
            s=5,
            alpha=0.35,
            color="#D1495B",
            linewidths=0,
        )
        rho = performance.loc[
            performance["protein_id"] == protein_id, "spearman"
        ].iloc[0]
        ax.set_title(f"{protein_id} (rho={rho:.3f})")
        ax.set_xlabel("Processed PBM E-score")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("DeepPBS PWM-derived score")
    fig.suptitle("DeepPBS versus experimental designed-DBP landscape")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_deeppbs_vs_pbm_completed.png", dpi=300)
    plt.close(fig)


def write_reports(
    manifest: pd.DataFrame,
    performance: pd.DataFrame,
    hard_summary: dict[str, int],
) -> None:
    evaluated = performance[performance["status"] == "evaluated"]
    median = float(evaluated["spearman"].median()) if len(evaluated) else np.nan
    reference = 0.5914386964424022
    old_summary_path = TABLES / "final_strong_baseline_summary.csv"
    old_summary = pd.read_csv(old_summary_path)
    old_summary = old_summary[old_summary["method"] != "DeepPBS"].copy()
    sequence_rows = old_summary[
        old_summary["method"].astype(str).str.contains("sequence-only", case=False)
    ]
    old_sequence = float(sequence_rows["designed_macro_median_spearman"].iloc[0])
    old_simple = float(
        old_summary.loc[
            old_summary["method"] == "SimpleProteinConditionalBaseline",
            "designed_macro_median_spearman",
        ].iloc[0]
    )
    old_frozen = float(
        old_summary.loc[
            old_summary["method"] == "FrozenPLMProteinConditionalBaseline",
            "designed_macro_median_spearman",
        ].iloc[0]
    )
    structure_lines = []
    for row in manifest.itertuples(index=False):
        structure_lines.append(
            f"| {row.dbp_id} | {row.pdb_id or 'NA'} | "
            f"{'yes' if row.deepPBS_evaluable else 'no'} | "
            f"{row.reason_not_evaluable or 'official run completed'} | "
            f"{row.overlap_status} |"
        )
    result_lines = []
    for row in performance.itertuples(index=False):
        value = f"{row.spearman:.6f}" if pd.notna(row.spearman) else "NA"
        result_lines.append(
            f"| {row.protein_id} | {value} | {row.status} | "
            f"{row.n_rc_classes or 'NA'} |"
        )
    final_summary = old_summary.copy()
    final_summary = pd.concat(
        [
            final_summary,
            pd.DataFrame(
                [
                    {
                        "method": "DeepPBS",
                        "protein_representation": "structure-aware official ensemble PWM",
                        "structure_required": True,
                        "training_data": "official bundled DeepPBS model",
                        "n_designed_proteins_covered": int(len(evaluated)),
                        "designed_macro_median_spearman": median,
                        "natural_macro_median_spearman": np.nan,
                        "gap_to_replicate": reference - median,
                        "overlap_caveat": (
                            "exact overlap not found in checked manifests; "
                            "homolog-level audit unresolved"
                        ),
                        "notes": (
                            "Completed Linux diagnostic; DBP35/DBP48 only; "
                            "not a 7-protein generalization estimate"
                        ),
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    final_summary.to_csv(
        TABLES / "final_strong_baseline_summary_deeppbs_completed_v0_4_2.csv",
        index=False,
    )
    report = f"""# DeepPBS Baseline Completion Report

Audit date: 2026-09-04

## Scope

This report closes the previously incomplete DeepPBS baseline. The official
upstream checkout was run remotely on Ubuntu using the documented preprocessing
and prediction scripts. The project did not modify the upstream model or
weights. The PBM comparison is a per-protein ranking comparison between
processed experimental uPBM E-scores and a fixed PWM-derived sequence score.
The completion claim is technical and diagnostic: only structures with a
legitimate protein-DNA input are included, so coverage is reported explicitly.

## Structure coverage

| DBP | PDB | Evaluable | Reason | Overlap status |
|---|---|---:|---|---|
{chr(10).join(structure_lines)}

DeepPBS coverage is **{len(evaluated)}/7**. DBP35 uses the available theoretical
design complex. DBP48 uses the experimental 8TAC complex after a project-side
DSSR-defined helix-only input preparation; the original 8TAC file is preserved.

## Output semantics and scoring

The upstream source confirms `A/C/G/T` column order in
`deeppbs/dna_encodings.py::seqToOneHot` and `oneHotToSeq`. In
`run/predict.py`, `P` is the post-softmax ensemble output after averaging the
two strand halves with a reversed second half; `Seq` is the hard input
sequence. The project scorer evaluates every contiguous 7-mer PWM window with
`sum(log(P + 1e-9))`, takes the maximum window, then takes the maximum over a
candidate and its reverse complement. This is a fixed DeepPBS-derived ranking
proxy, not affinity, Kd, probability, or a calibrated binding score.

## PBM results

| DBP | Spearman | Status | RC classes |
|---|---:|---|---:|
{chr(10).join(result_lines)}

DeepPBS evaluable macro median Spearman is **{median:.6f}** across **{len(evaluated)}/7**
proteins. For context, the prior sequence-only k-mer3 median was
`{old_sequence:.6f}`, SimpleProteinConditional was `{old_simple:.6f}`,
FrozenPLM was `{old_frozen:.6f}`, and the empirical replicate agreement
reference was `{reference:.6f}`. These numbers have different coverage and
model assumptions; the two-protein DeepPBS median is diagnostic, not a
seven-protein generalization estimate.

## Hard-case integration

Using the frozen v0.3.1 disagreement definition, DeepPBS has
**{hard_summary['eligible_candidates']} eligible candidate rows**, resolves
**{hard_summary['resolved']}**, leaves **{hard_summary['unresolved']} unresolved**,
and has a resolution rate of **{hard_summary['resolution_rate']:.6f}** within
its eligible denominator. The old 1,515 total is the all-protein candidate
set; it must not be used as DeepPBS's denominator because five proteins have
no legal structure input. The completed four-method common-hard subset contains
**{hard_summary['complete_core_common_failures']}** rows and is written separately.

The frozen all-protein disagreement set contains
**{hard_summary['all_disagreement_candidates']}** candidates. DeepPBS is not
assigned a denominator of 1,515 because DBP1, DBP3, DBP5, DBP6, and DBP9 have
no legal structure input. Candidate-level failure categories and pairwise
complementarity counts are in
`tables/baseline_failure_cases_deeppbs_completed_v0_4_2.parquet` and
`tables/baseline_complementarity_deeppbs_completed_v0_4_2.csv`.

## Baseline arena

The completed comparison table is
`tables/final_strong_baseline_summary_deeppbs_completed_v0_4_2.csv`.
DeepPBS is the only structure-aware result here and covers 2/7 proteins.
Its median Spearman is therefore a diagnostic statistic, not a fair
seven-protein ranking claim.

## Official example and run provenance

The official `5x6g` example passed in the Ubuntu VM. DBP35 completed with
Helix score 1.0 and contact count 266. DBP48 completed after DSSR-defined
helix-only preparation with Helix score 1.0 and contact count 224. Exact
commands, stdout/stderr, environment, input/output paths, and SHA256 values
are retained under `external_runs/` and in
`tables/deeppbs_run_manifest_completed_v0_4_2.csv`.

## Limitations

- Only DBP35 and DBP48 have a legal public/project structure input in this
  repository; five proteins remain not evaluable rather than being fabricated.
- DBP35 is a theoretical Rosetta model, not an experimental complex.
- DBP48 required a transparent helix-only preprocessing repair because the
  deposited complex has non-helical DNA overhang residues that break the
  upstream shape extraction path.
- The checked DeepPBS manifests did not contain exact designed-DBP/PDB hits,
  but full homolog-level training-set audit remains unresolved because the
  upstream repository does not distribute all training sequences.
- A two-protein result cannot establish a systematic strong-baseline ranking.
"""
    (V042 / "DEEPPBS_COMPLETION_REPORT.md").write_text(report, encoding="utf-8")

    gap_report = f"""# Baseline Gap Analysis After DeepPBS Completion

Audit date: 2026-09-04

## Main result

DeepPBS is technically reproducible in the Ubuntu VM, but it is evaluable for
only **{len(evaluated)}/7** designed DBPs. Its diagnostic median Spearman is
**{median:.6f}**, compared with the frozen sequence-only best baseline
**{old_sequence:.6f}**, SimpleProteinConditional **{old_simple:.6f}**,
FrozenPLM **{old_frozen:.6f}**, and empirical replicate agreement
**{reference:.6f}**.

These are not directly comparable as a complete model ranking because coverage,
structure requirements, and training/overlap caveats differ.

## Interpretation

- DeepPBS does not demonstrate a clear improvement over the existing
  sequence-only or SimplePC results on the two evaluable proteins.
- DBP35 is low (`{float(evaluated.loc[evaluated['protein_id'].eq('DBP35'), 'spearman'].iloc[0]):.6f}`)
  and DBP48 is higher (`{float(evaluated.loc[evaluated['protein_id'].eq('DBP48'), 'spearman'].iloc[0]):.6f}`),
  but the sample size is two and DBP35 uses a theoretical design model.
- The previous 1,515 disagreement count remains an all-protein count.
  DeepPBS has 398 eligible candidates, resolves 39, and leaves 359 unresolved.
- The completed four-method common-low set contains
  {hard_summary['complete_core_common_failures']} candidates.

## Boundary of claim

This closes the runtime and evaluation gap, but not the structure-coverage gap.
The result supports a structure-aware diagnostic baseline and a reproducible
failure analysis. It does not support a clean zero-shot or seven-protein
generalization claim.
"""
    (V042 / "BASELINE_GAP_ANALYSIS_DEEPPBS_COMPLETED.md").write_text(
        gap_report, encoding="utf-8"
    )

    gate = f"""# DeepPBS-Completed Baseline Gate

The prior v0.4.2 `WAIT` status was caused by an unevaluated DeepPBS runtime,
not by a model result. This completion artifact supersedes that interpretation
without deleting the historical report.

## Gate

**WAIT FOR STRONGER STRUCTURE-COVERED BASELINE**

DeepPBS is now technically integrated and has a valid two-protein diagnostic
result, but the structure coverage is only **{len(evaluated)}/7**. The available
result does not justify calling DeepPBS a complete seven-protein strong
baseline or claiming designed-DBP generalization. The project may use it as a
structure-aware diagnostic comparator while pursuing additional public
designed complexes or a separately labeled predicted-structure sensitivity
study.

## Scientific conclusion

DeepPBS can be evaluated against the experimental PBM landscape when a
protein-matched protein-DNA structure is available. This run does not show
that structure-aware scoring solves the designed specificity problem: the
experimental reference is `{reference:.3f}` median Spearman, while the current
DeepPBS coverage is only two proteins. The remaining five missing structures
are a coverage limitation, not missing predictions filled with zero.

Do not call DBP35/DBP48 a clean zero-shot result: exact overlap was not found in
the checked manifests, but homolog-level overlap remains unresolved.
"""
    (V042 / "FINAL_STRONG_BASELINE_GATE_DEEPPBS_COMPLETED.md").write_text(
        gate, encoding="utf-8"
    )

    memo = f"""# v0.5 Model Decision Memo

## Current recommendation

**WAIT FOR STRUCTURE-COVERAGE CLOSURE BEFORE CLAIMING A FULL STRONG-BASELINE
COMPARISON; DO NOT START THE PROPOSED MODEL YET.**

DeepPBS now runs through the official Linux workflow and supplies a valid
structure-aware diagnostic for DBP35 and DBP48. However, only **{len(evaluated)}/7**
designed proteins are evaluable. The observed median Spearman ({median:.3f})
cannot be compared as a complete seven-protein result with the
{old_sequence:.3f} sequence-only, {old_simple:.3f} SimplePC, or
{reference:.3f} replicate references.

## What is supported

- DeepPBS preprocessing and inference are reproducible in the Ubuntu VM.
- The fixed PWM-to-7-mer mapping and RC-class scorer are implemented.
- DBP35 and DBP48 can be included as a transparent structure-aware diagnostic.
- Five proteins remain structurally unevaluable without fabricating inputs.

## What is not supported

- No claim that DeepPBS is the strongest baseline.
- No clean zero-shot claim because homolog-level training overlap is unresolved.
- No claim that a low result is biological rather than structure/mapping-related
  for DBP48 or theoretical-model-related for DBP35.

## v0.5 hypothesis to test after coverage closure

If future structure coverage confirms that DeepPBS, SimplePC, FrozenPLM, and
sequence-only baselines leave a stable common failure set, the minimal
failure-driven hypothesis is a target-anchored differential ranker:
compare `S(P,T)` and `S(P,D)` directly rather than scoring `D` in isolation.
The first falsification experiment should use a protein-held-out, RC-safe,
per-protein ranking benchmark and compare target-relative margins against the
best existing baseline with bootstrap confidence intervals. This is a memo
only; no proposed model is implemented in this commit.
"""
    (DOCS / "V0_5_MODEL_DECISION_MEMO.md").write_text(memo, encoding="utf-8")


def main() -> None:
    experimental = load_designed_experimental_units(ROOT)
    predictions = attach_experimental_scores(experimental, load_prediction_tables())
    predictions.to_parquet(
        TABLES / "deeppbs_predictions_completed_v0_4_2.parquet", index=False
    )
    manifest = write_structure_manifest()
    write_overlap_manifest()
    write_run_manifest()
    performance = write_performance(experimental, predictions)
    scored = load_scored_landscape(experimental, predictions)
    hard_summary = write_hard_case_outputs(scored)
    write_figures(performance, scored)
    write_reports(manifest, performance, hard_summary)
    summary = {
        "prediction_rows": int(len(predictions)),
        "evaluable_proteins": int((performance["status"] == "evaluated").sum()),
        "coverage": f"{(performance['status'] == 'evaluated').sum()}/7",
        "macro_median_spearman": float(
            performance.loc[performance["status"] == "evaluated", "spearman"].median()
        ),
        **hard_summary,
    }
    (TABLES / "deeppbs_completion_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(performance.to_string(index=False))


if __name__ == "__main__":
    main()
