from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_v0_4_2_construct_audit_marks_all_natural_proteins_unknown():
    audit = pd.read_csv(ROOT / "metadata" / "v0_4_2" / "natural_pbm_construct_audit.csv")
    assert audit["protein_id"].nunique() == 57
    assert audit["experimental_construct_known"].eq(False).all()
    assert audit["use_for_primary_training"].eq(False).all()


def test_v0_4_2_assay_aligned_benchmark_is_empty():
    aligned = pd.read_parquet(ROOT / "data" / "processed" / "v0_4_2" / "natural_pbm_assay_aligned_v0_4_2.parquet")
    assert aligned.shape[0] == 0


def test_v0_4_2_deeppbs_official_example_not_run_on_host():
    stderr = (ROOT / "results" / "v0_4_2" / "external_runs" / "deeppbs_official_example" / "stderr.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    assert "OFFICIAL_EXAMPLE_NOT_RUN_HOST_NO_LINUX_RUNTIME" in stderr


def test_v0_4_2_frozen_plm_schema_and_scores():
    perf = pd.read_csv(ROOT / "results" / "v0_4_2" / "tables" / "frozen_plm_performance.csv")
    macro = pd.read_csv(ROOT / "results" / "v0_4_2" / "tables" / "frozen_plm_performance_macro.csv")
    assert set(perf["dataset"]) == {"natural_test", "designed_external"}
    assert perf["protein_id"].nunique() == 16
    assert macro[macro["metric"].eq("spearman") & macro["dataset"].eq("designed_external")]["median"].iloc[0] > 0.0
    assert macro[macro["metric"].eq("spearman") & macro["dataset"].eq("natural_test")]["median"].iloc[0] > 0.0


def test_v0_4_2_summary_records_wait_gate():
    gate = (ROOT / "results" / "v0_4_2" / "FINAL_STRONG_BASELINE_GATE.md").read_text(encoding="utf-8")
    assert "WAIT - BENCHMARK STILL INCOMPLETE" in gate


def test_v0_4_2_no_absolute_windows_paths_in_new_manifests():
    for rel in [
        "metadata/v0_4_2/deeppbs_weight_manifest_v2.csv",
        "metadata/v0_4_2/designed_structure_manifest_v2.csv",
        "metadata/v0_4_2/natural_pbm_construct_audit.csv",
    ]:
        text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
        assert "E:\\" not in text
        assert "C:\\" not in text

