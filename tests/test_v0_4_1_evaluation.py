from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_simple_pc_performance_is_per_protein_for_both_datasets():
    perf = pd.read_csv(ROOT / "results" / "v0_4_1" / "tables" / "simple_pc_performance.csv")
    assert set(perf["dataset"]) == {"natural_test", "designed_external"}
    assert perf[perf["dataset"] == "natural_test"]["protein_id"].nunique() >= 5
    assert perf[perf["dataset"] == "designed_external"]["protein_id"].nunique() == 7
    assert perf["spearman"].notna().all()


def test_simple_pc_macro_summary_keeps_dataset_dimension():
    macro = pd.read_csv(ROOT / "results" / "v0_4_1" / "tables" / "simple_pc_performance_macro.csv")
    assert {"dataset", "baseline", "metric", "median", "n_proteins"}.issubset(macro.columns)
    assert set(macro["dataset"]) == {"natural_test", "designed_external"}
    assert not macro[(macro["metric"] == "spearman") & (macro["dataset"] == "designed_external")].empty


def test_disagreement_resolution_uses_full_v0_3_1_candidate_counts():
    resolution = pd.read_csv(ROOT / "results" / "v0_4_1" / "tables" / "simple_pc_disagreement_resolution.csv")
    assert int(resolution["n_v0_3_1_disagreement_candidates"].sum()) == 1515
    assert resolution["n_evaluable_candidates"].eq(resolution["n_v0_3_1_disagreement_candidates"]).all()
    assert resolution["n_resolved"].le(resolution["n_evaluable_candidates"]).all()


def test_final_gate_does_not_claim_deeppbs_success():
    text = (ROOT / "results" / "v0_4_1" / "FINAL_MODEL_DEVELOPMENT_GATE.md").read_text(encoding="utf-8")
    assert "WAIT FOR STRONGER BASELINE" in text
    assert "DeepPBS was not fairly run" in text
