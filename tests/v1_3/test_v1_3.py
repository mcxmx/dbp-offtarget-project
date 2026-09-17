from pathlib import Path

import numpy as np
import pandas as pd

import scripts.v1_3.run_v1_3 as v13


ROOT = Path(__file__).resolve().parents[2]


def test_competition_tidy_table_schema_and_counts() -> None:
    table = pd.read_parquet(ROOT / "data/processed/v1_3_competition_mutations.parquet")
    required = {
        "protein_id", "assay", "target_sequence", "position", "wt_base",
        "mutant_base", "replicate", "raw_measurement", "normalized_effect",
        "source_file", "development_status", "holdout_status",
    }
    assert required <= set(table.columns)
    assert len(table) == 414
    assert table.protein_id.nunique() == 10
    assert table.groupby(["protein_id", "position"]).size().eq(3).all()
    assert not table.duplicated(
        ["protein_id", "position", "wt_base", "mutant_base"]
    ).any()
    np.testing.assert_allclose(table.normalized_effect, -table.raw_measurement)


def test_holdout_status_was_frozen_before_method_results() -> None:
    table = pd.read_parquet(ROOT / "data/processed/v1_3_competition_mutations.parquet")
    locked = set(table.loc[table.holdout_status == "LOCKED_HOLDOUT", "protein_id"])
    exposed = set(table.loc[table.holdout_status == "DEVELOPMENT_EXPOSED", "protein_id"])
    assert locked == {"DBP023", "DBP056", "DBP062"}
    assert exposed == set(v13.OLD_PROTEINS)
    assert locked.isdisjoint(set(pd.read_csv(
        ROOT / "results/v1_3/decomposition_per_protein.tsv", sep="\t"
    ).protein_id))


def test_primary_six_and_dbp048_sensitivity_only() -> None:
    assert v13.PRIMARY_PROTEINS == [
        "DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP035"
    ]
    per = pd.read_csv(ROOT / "results/v1_3/decomposition_per_protein.tsv", sep="\t")
    assert set(per.loc[per.analysis_status == "PRIMARY", "protein_id"]) == set(v13.PRIMARY_PROTEINS)
    assert set(per.loc[per.analysis_status == "SENSITIVITY_ONLY", "protein_id"]) == {"DBP048"}
    summary = pd.read_csv(ROOT / "results/v1_3/decomposition_summary.tsv", sep="\t")
    assert set(summary.analysis_set) == {"primary_six", "all_seven_descriptive"}
    assert set(summary.loc[summary.analysis_set == "primary_six", "n_proteins"]) == {6}


def test_decomposition_coverage_and_figure2_traceability() -> None:
    per = pd.read_csv(ROOT / "results/v1_3/decomposition_per_protein.tsv", sep="\t")
    fig = pd.read_csv(ROOT / "results/v1_3/figure2_data.tsv", sep="\t")
    assert set(per.method) == {"DeepPBS", "NA-MPNN", "PBM_experimental"}
    assert per.groupby("method").size().to_dict() == {
        "DeepPBS": 7, "NA-MPNN": 7, "PBM_experimental": 7
    }
    assert set(per.loc[per.method == "DeepPBS", "n_mutations"]) == {39, 42}
    assert set(per.loc[per.method != "DeepPBS", "n_mutations"]) == {21}
    assert {
        "position_mean_signed_spearman", "position_max_abs_spearman",
        "within_kendall_style_concordance",
    } <= set(per.columns)
    np.testing.assert_allclose(
        per.within_kendall_style_concordance,
        2 * per.within_pairwise_accuracy - 1,
        atol=1e-15,
    )
    keys = ["method", "protein_id"]
    merged = fig.merge(per, on=keys, suffixes=("_fig", "_per"), validate="one_to_one")
    assert len(merged) == len(per) == 21
    np.testing.assert_allclose(merged.global_spearman_fig, merged.global_spearman_per)
    np.testing.assert_allclose(merged.position_spearman_fig, merged.position_spearman_per)
    np.testing.assert_allclose(
        merged.within_position_residual_spearman,
        merged.within_residual_spearman,
    )
    assert set(fig.global_chance) == {0.0}
    assert set(fig.position_chance) == {0.0}
    assert set(fig.within_pairwise_chance) == {0.5}


def test_permutation_null_outputs_and_loo() -> None:
    nulls = pd.read_csv(ROOT / "results/v1_3/permutation_nulls.tsv", sep="\t")
    distribution = pd.read_csv(ROOT / "results/v1_3/permutation_distributions.tsv", sep="\t")
    assert len(nulls) == 3 * 7 * 3
    assert len(distribution) == 3 * 7 * 3 * v13.N_PERM
    assert set(nulls.metric) == {
        "position_spearman", "within_residual_spearman", "within_pairwise_accuracy"
    }
    assert nulls.empirical_p.between(0, 1, inclusive="both").all()
    assert set(nulls.n_permutations) == {v13.N_PERM}
    loo = pd.read_csv(ROOT / "results/v1_3/leave_one_protein_out.tsv", sep="\t")
    assert len(loo) == 3 * 6 * 4
    assert set(loo.n_proteins) == {5}


def test_permutation_implementation_is_deterministic(monkeypatch, tmp_path: Path) -> None:
    comp = pd.read_parquet(ROOT / "data/processed/v1_3_competition_mutations.parquet")
    comp = comp[comp.protein_id == "DBP001"]
    monkeypatch.setattr(v13, "OUT", tmp_path)
    monkeypatch.setattr(v13, "N_PERM", 16)
    first = v13.permutation_nulls(comp, pd.DataFrame())
    first_distribution = pd.read_csv(tmp_path / "permutation_distributions.tsv", sep="\t")
    second = v13.permutation_nulls(comp, pd.DataFrame())
    second_distribution = pd.read_csv(tmp_path / "permutation_distributions.tsv", sep="\t")
    pd.testing.assert_frame_equal(first, second)
    pd.testing.assert_frame_equal(first_distribution, second_distribution)


def test_contract_records_unavailable_dbp35opt() -> None:
    contract = (ROOT / "V1_3_EVALUATION_CONTRACT.md").read_text(encoding="utf-8")
    assert "FROZEN 2026-09-17" in contract
    assert "DBP048 remains `SENSITIVITY_ONLY`" in contract
    assert "DBP35opt" in contract and "NOT EVALUABLE" in contract


def test_historical_shuffle_caches_are_reused_and_standardized() -> None:
    detail = pd.read_csv(ROOT / "results/v1_3/cached_shuffle_diagnostics.tsv", sep="\t")
    summary = pd.read_csv(ROOT / "results/v1_3/cached_shuffle_summary.tsv", sep="\t")
    assert set(detail.shuffle_type) == {"protein", "target"}
    assert set(detail.model) == {"M3"}
    assert set(detail.development_status) == {"development_exposed_historical_v0_6"}
    assert set(summary.model) == {"M0", "M1", "M1c", "M2", "M3"}
    assert "target_conditioning_effect_size_median" in summary.columns
