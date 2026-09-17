from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "v0_7_benchmark"


def test_required_v07_artifacts_exist():
    required = [
        ROOT / "docs/v0_7_benchmark/BENCHMARK_PROTOCOL.md",
        ROOT / "docs/v0_7_benchmark/HARD_CASE_ANALYSIS.md",
        ROOT / "docs/v0_7_benchmark/NATURAL_DESIGNED_DOMAIN_SHIFT.md",
        ROOT / "docs/v0_7_benchmark/V0_7_SCIENTIFIC_ASSESSMENT.md",
        OUT / "benchmark_master.csv",
        OUT / "performance_conditionality.csv",
        OUT / "hard_case_taxonomy.csv",
        OUT / "protein_sample_size_simulation.csv",
        OUT / "natural_designed_shift.csv",
        ROOT / "metadata/v0_7_external/candidate_datasets.csv",
        ROOT / "metadata/v0_7_benchmark/sota_method_manifest.csv",
    ]
    assert all(path.exists() for path in required)


def test_master_table_preserves_protocol_and_exposure_fields():
    frame = pd.read_csv(OUT / "benchmark_master.csv")
    required = {
        "dataset", "training_proteins", "held_out_proteins", "assay", "dna_kmer_length",
        "split_protocol", "exposure_status", "model_selection_exposure", "median_spearman",
        "per_protein_spearman", "random_seeds", "protein_conditionality_pass",
        "target_conditionality_pass",
    }
    assert required.issubset(frame.columns)
    assert frame["exposure_status"].astype(str).str.contains("exposed", case=False).all()
    assert not frame["comparability_status"].eq("independent external validation").any()
    assert frame["median_spearman"].notna().all()


def test_conditionality_axis_records_collapse_thresholds():
    frame = pd.read_csv(OUT / "performance_conditionality.csv")
    assert {"stage", "method", "median_heldout_spearman", "protein_shuffle_correlation", "target_shuffle_correlation", "four_state"}.issubset(frame.columns)
    m3 = frame.loc[(frame.stage == "v0.6") & frame.method.eq("M3")].iloc[0]
    assert m3.protein_shuffle_correlation > 0.995
    assert m3.target_shuffle_correlation > 0.995
    assert m3.four_state == "low_performance+collapsed_conditioning"


def test_hard_case_taxonomy_contains_registered_and_descriptive_strata():
    frame = pd.read_csv(OUT / "hard_case_taxonomy.csv")
    categories = set(frame.category)
    assert {"all_model_failure", "sequence_only_rescued", "conditional_m1c_rescued", "joint_conditional_rescued"}.issubset(categories)
    assert {"motif_near_target", "highly_dissimilar_off_target", "high_experimental_low_predicted", "low_experimental_high_predicted"}.issubset(categories)
    assert int(frame.loc[frame.category.eq("all_model_failure"), "count"].iloc[0]) == 871
    assert {"fraction_of_reference_1515", "fraction_of_denominator", "denominator_n"}.issubset(frame.columns)
    near = frame.loc[frame.category.eq("motif_near_target")].iloc[0]
    assert near.denominator_n == 871
    assert abs(near.fraction_of_denominator - near["count"] / near.denominator_n) < 1e-12


def test_simulation_is_explicitly_non_empirical():
    frame = pd.read_csv(OUT / "protein_sample_size_simulation.csv")
    assert frame.simulation_only.eq(True).all()
    assert set(frame.n_training_proteins) == {2, 3, 4, 5, 10, 20, 40, 80}
    assert frame.n_dna_per_protein.nunique() == 1


def test_external_firewall_has_no_quantitative_labels_and_candidate_status():
    frame = pd.read_csv(ROOT / "metadata/v0_7_external/candidate_datasets.csv")
    assert frame.exposure_status.str.contains("UNEXPOSED_CONFIRMATORY_CANDIDATE").any()
    assert "quantitative_labels" not in frame.columns
    assert not frame.study.str.contains("specificity result", case=False).any()


def test_figure_five_is_blocked_until_independent_cohort():
    frame = pd.read_csv(OUT / "figure_data_interface_status.csv")
    row = frame.loc[frame.figure.eq("Figure 5")].iloc[0]
    assert row.status == "BLOCKED"
