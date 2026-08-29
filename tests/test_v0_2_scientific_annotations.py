import math
from pathlib import Path

import pandas as pd

from src.models.base_binding_model import SequenceProxyBaseline


ROOT = Path(__file__).resolve().parents[1]


def as_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def test_pdb_structure_does_not_imply_quantitative_specificity_ground_truth():
    curated = pd.read_csv(ROOT / "data" / "processed" / "dbp_target_pairs_all_curated.csv")
    pdb_rows = curated[curated["source_type"] == "PDB"]
    assert not pdb_rows.empty
    assert as_bool_series(pdb_rows["has_structural_cognate"]).all()
    assert not as_bool_series(pdb_rows["has_quantitative_specificity_ground_truth"]).any()


def test_v0_2_single_mutant_counts_match_theory():
    pairs = pd.read_csv(ROOT / "data" / "processed" / "dbp_target_pairs_v0_2.csv")
    singles = pd.read_csv(ROOT / "data" / "processed" / "single_mutants_v0_2.csv")
    observed = singles.groupby("pair_id").size().to_dict()
    for _, pair in pairs.iterrows():
        assert observed[pair["pair_id"]] == 3 * int(pair["dna_length"])


def test_v0_2_double_mutant_counts_match_theory():
    pairs = pd.read_csv(ROOT / "data" / "processed" / "dbp_target_pairs_v0_2.csv")
    doubles = pd.read_csv(ROOT / "data" / "processed" / "double_mutants_v0_2.csv")
    observed = doubles.groupby("pair_id").size().to_dict()
    for _, pair in pairs.iterrows():
        expected = 9 * math.comb(int(pair["dna_length"]), 2)
        assert observed[pair["pair_id"]] == expected


def test_sequence_proxy_baseline_is_not_protein_conditioned():
    model = SequenceProxyBaseline()
    assert model.is_protein_conditioned is False


def test_core_benchmark_excludes_non_specific_and_guide_dependent():
    pairs = pd.read_csv(ROOT / "data" / "processed" / "dbp_target_pairs_v0_2.csv")
    forbidden = {"non_specific", "guide_dependent"}
    assert not (set(pairs["sequence_specificity_class"]) & forbidden)


def test_curated_provenance_has_source_id_or_url():
    curated = pd.read_csv(ROOT / "data" / "processed" / "dbp_target_pairs_all_curated.csv")
    has_id_or_url = curated["source_id"].fillna("").astype(str).ne("") | curated["source_url"].fillna("").astype(str).ne("")
    assert has_id_or_url.all()

