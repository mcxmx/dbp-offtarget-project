from pathlib import Path

import numpy as np
import pandas as pd

from scripts.v1_1.run_v1_1 import PROTEINS, canonical_landscape_from_pwm, decode_one_hot


ROOT = Path(__file__).resolve().parents[2]


def test_canonical_landscape_is_8192_rows():
    P = np.full((7, 4), 0.25, dtype=float)
    out = canonical_landscape_from_pwm(P, list(range(7)))
    assert len(out) == 8192
    assert out.canonical_7mer.is_unique


def test_deeppbs_bundle_has_all_seven_valid_npz():
    root = ROOT / "data/raw/v1_1_deeppbs/frozen_design7_20260914/predictions"
    for protein in PROTEINS:
        data = np.load(root / f"{protein}.npz_predict.npz", allow_pickle=False)
        assert set(data.files) == {"P", "Seq"}
        assert data["P"].shape == data["Seq"].shape
        assert data["P"].shape[1] == 4
        assert np.isfinite(data["P"]).all()
        assert np.allclose(data["P"].sum(axis=1), 1, atol=2e-5)


def test_competition_source_is_parsed_without_fabricated_replicates():
    status = pd.read_csv(ROOT / "results/v1_1/04_competition/competition_source_status.tsv", sep="\t")
    assert status.iloc[0].status == "PARSED_MOESM16"
    mutations = pd.read_csv(ROOT / "results/v1_1/04_competition/competition_mutations.tsv", sep="\t")
    assert len(mutations) == 288
    assert set(mutations.protein) == set(PROTEINS)
    assert not mutations.replicate_level_available.any()


def test_competition_has_three_substitutions_per_position():
    mutations = pd.read_csv(ROOT / "results/v1_1/04_competition/competition_mutations.tsv", sep="\t")
    counts = mutations.groupby(["protein", "assay_position"]).mut_base.nunique()
    assert counts.eq(3).all()
    assert mutations.experimental_effect.equals(-mutations.normalized_value)


def test_dbp048_uses_sequence_c_and_remains_sensitivity_only():
    mapping = pd.read_csv(ROOT / "results/v1_1/03_register_mapping/register_summary.tsv", sep="\t")
    row = mapping[mapping.protein == "DBP048"].iloc[0]
    assert row.assay_sequence == "CGACACCTGACGCG"
    assert row.mapping_status == "AMBIGUOUS"


def test_competition_summary_has_no_missing_source_status():
    summary = (ROOT / "results/v1_1/v1_1_summary.md").read_text(encoding="utf-8")
    assert "MISSING_LOCAL_SOURCE" not in summary
    assert "ASSAY-ALIGNED DIAGNOSTIC COMPLETE" in summary


def test_global_projection_has_full_per_protein_universe():
    df = pd.read_csv(ROOT / "results/v1_1/09_global_pwm_projection/unified_global_landscapes.tsv", sep="\t")
    assert set(df.protein) == set(PROTEINS)
    assert set(df.groupby(["method", "protein"]).size()) == {8192}


def test_no_v1_0_outputs_are_written_by_v1_1():
    assert not list((ROOT / "results/v1_0_structure").glob("*v1_1*"))
