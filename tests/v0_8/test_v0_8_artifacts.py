from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def test_required_v08_artifacts_exist():
    required = [
        "docs/v0_8_sota/SOTA_BENCHMARK_PROTOCOL.md",
        "docs/v0_8_sota/AF3_FEASIBILITY_AUDIT.md",
        "docs/v0_8_sota/V0_8_SOTA_RESULTS.md",
        "docs/v0_8_external/EXTERNAL_COHORT_SELECTION_PROTOCOL.md",
        "docs/v0_8_external/PREREGISTERED_VALIDATION_PLAN.md",
        "docs/paper/PAPER_OUTLINE_V0_8.md",
        "metadata/v0_8_sota/method_execution_manifest.csv",
        "metadata/v0_8_external/candidate_datasets.csv",
        "metadata/v0_8_external/FROZEN_VALIDATION_MANIFEST.yaml",
        "results/v0_8_sota/sota_benchmark_master.csv",
        "results/v0_8_sota/experimental_reference.csv",
    ]
    assert all((ROOT / path).exists() for path in required)


def test_structure_metrics_are_coverage_aware():
    master = pd.read_csv(ROOT / "results/v0_8_sota/sota_benchmark_master.csv")
    for method in ("DeepPBS", "NA-MPNN"):
        row = master.loc[master.method == method].iloc[0]
        assert row.coverage == "2/7"
        assert row.protein_unit_n == 2
        assert "partial coverage" in row.comparability_status


def test_sequence_baseline_is_one_metric_per_protein():
    master = pd.read_csv(ROOT / "results/v0_8_sota/sota_benchmark_master.csv")
    row = master.loc[master.method == "sequence_kmer3"].iloc[0]
    assert row.protein_unit_n == 7
    assert abs(row.median_spearman - 0.2320822289869546) < 1e-12
    assert row.per_protein_spearman.count("DBP") == 7


def test_external_firewall_is_label_locked():
    candidates = pd.read_csv(ROOT / "metadata/v0_8_external/candidate_datasets.csv")
    assert candidates.label_access.eq(False).all()
    assert "NO SUITABLE PUBLIC CONFIRMATORY COHORT" in (
        ROOT / "docs/v0_8_external/EXTERNAL_COHORT_SELECTION_PROTOCOL.md"
    ).read_text(encoding="utf-8")
    manifest = (ROOT / "metadata/v0_8_external/FROZEN_VALIDATION_MANIFEST.yaml").read_text(encoding="utf-8")
    assert "label_access: false" in manifest


def test_af3_is_explicitly_not_run():
    metrics = pd.read_csv(ROOT / "results/v0_8_sota/af3_contact_metrics.csv")
    assert metrics.empty
    assert "FEASIBILITY_ONLY" in (ROOT / "docs/v0_8_sota/AF3_FEASIBILITY_AUDIT.md").read_text(encoding="utf-8")

