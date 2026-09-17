from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

import scripts.v1_2.run_v1_2 as v12
from src.sequence_equivalence import canonical_rc, reverse_complement


ROOT = Path(__file__).resolve().parents[2]
EXPECTED_V11 = {
    "DBP001": 0.2271255060728745,
    "DBP003": 0.18238866396761136,
    "DBP005": 0.34916133214488293,
    "DBP006": 0.7447532614861032,
    "DBP009": 0.5655133295519001,
    "DBP035": 0.3583988331577668,
    "DBP048": 0.049833887043189376,
}


def _deep() -> pd.DataFrame:
    return pd.read_csv(
        ROOT / "results/v1_1/05_assay_aligned_eval/deeppbs_per_mutation_predictions.tsv",
        sep="\t",
    )


def _pbm() -> pd.DataFrame:
    return pd.read_csv(
        ROOT / "results/v1_1/06_assay_agreement/pbm_competition_per_mutation.tsv",
        sep="\t",
    )


def _prepare_out(monkeypatch, tmp_path: Path, proteins: list[str]) -> None:
    monkeypatch.setattr(v12, "OUT", tmp_path)
    monkeypatch.setattr(v12, "PROTEINS", proteins)
    for directory in ["02_permutation_nulls", "07_assay_triangle", "10_position_robustness"]:
        (tmp_path / directory).mkdir(parents=True, exist_ok=True)


def test_exact_v11_deeppbs_reproduction() -> None:
    deep = _deep()
    assert set(deep.protein) == set(EXPECTED_V11)
    for protein, expected in EXPECTED_V11.items():
        group = deep[deep.protein == protein]
        observed = v12.safe_spearman(group.predicted_delta_logP, group.experimental_effect)
        assert np.isclose(observed, expected, atol=1e-12)


def test_permutation_nulls_are_deterministic(monkeypatch, tmp_path: Path) -> None:
    _prepare_out(monkeypatch, tmp_path, ["DBP001"])
    monkeypatch.setattr(v12, "N_PERMUTATIONS", 64)
    frame = _deep().query("protein == 'DBP001'")

    summary_a, distribution_a = v12.permutation_nulls(frame)
    summary_b, distribution_b = v12.permutation_nulls(frame)

    pd.testing.assert_frame_equal(summary_a, summary_b)
    pd.testing.assert_frame_equal(distribution_a, distribution_b)


def test_permutation_value_helper_is_deterministic() -> None:
    group = _deep().query("protein == 'DBP001'")
    first = v12.generate_permutation_values(group, "position_preserving_shuffle", 64, 2718)
    second = v12.generate_permutation_values(group, "position_preserving_shuffle", 64, 2718)
    np.testing.assert_array_equal(first, second)


def test_position_preserving_shuffle_does_not_cross_positions(
    monkeypatch, tmp_path: Path
) -> None:
    _prepare_out(monkeypatch, tmp_path, ["DBP001"])
    monkeypatch.setattr(v12, "N_PERMUTATIONS", 128)
    # Each position has one fixed prediction value. A within-position shuffle
    # must therefore leave the vector and its correlation unchanged exactly.
    frame = pd.DataFrame(
        {
            "protein": ["DBP001"] * 9,
            "assay_position": np.repeat([1, 2, 3], 3),
            "predicted_delta_logP": np.repeat([-2.0, 0.0, 3.0], 3),
            "experimental_effect": [-3, -2, -1, -1, 0, 1, 1, 2, 3],
        }
    )

    summary, distribution = v12.permutation_nulls(frame)
    observed = summary.loc[
        summary.null_type == "position_preserving_shuffle", "observed_spearman"
    ].iloc[0]
    within = distribution.loc[
        distribution.null_type == "position_preserving_shuffle", "spearman"
    ]
    assert np.allclose(within, observed)
    assert distribution.loc[
        distribution.null_type == "all_mutation_shuffle", "spearman"
    ].nunique() > 1


def test_triangle_uses_exact_matched_mutation_keys(monkeypatch, tmp_path: Path) -> None:
    _prepare_out(monkeypatch, tmp_path, list(v12.PROTEINS))
    deep, pbm = _deep(), _pbm()
    triangle, matched = v12.assay_triangle(deep, pbm)
    keys = ["protein", "assay_position", "wt_base", "mut_base"]

    assert not matched.duplicated(keys).any()
    assert (triangle.n_matched == 21).all()
    expected = deep[keys].merge(pbm[keys], on=keys, how="inner")
    pd.testing.assert_frame_equal(
        matched[keys].sort_values(keys).reset_index(drop=True),
        expected[keys].sort_values(keys).reset_index(drop=True),
    )


def test_dbp048_is_excluded_from_primary_six() -> None:
    assert v12.PRIMARY == ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP035"]
    assert "DBP048" not in v12.PRIMARY
    summary = pd.read_csv(ROOT / "results/v1_2/10_position_robustness/summary.tsv", sep="\t")
    assert summary.loc[summary.protein == "DBP048", "status"].iloc[0] == "SENSITIVITY_ONLY"


def test_constant_vector_correlations_are_explicitly_undefined() -> None:
    assert np.isnan(v12.safe_spearman([1, 1, 1], [1, 2, 3]))
    assert np.isnan(v12.safe_spearman([1, 2, 3], [4, 4, 4]))
    assert np.isnan(v12.safe_pearson([1, 1, 1], [1, 2, 3]))


def test_leave_one_position_out_removes_exactly_three_mutations(
    monkeypatch, tmp_path: Path
) -> None:
    _prepare_out(monkeypatch, tmp_path, ["DBP001"])
    group = _deep().query("protein == 'DBP001'")
    detail, summary = v12.position_robustness(group)

    assert len(group) == 39
    assert len(detail) == 13
    assert set(detail.n_remaining) == {36}
    assert set(detail.dropped_position) == set(group.assay_position)
    assert summary.n_positions.iloc[0] == 13


def test_rc_canonicalization_retains_8192_equivalence_classes() -> None:
    kmers = ["".join(chars) for chars in product("ACGT", repeat=7)]
    canonical = {canonical_rc(kmer) for kmer in kmers}
    assert len(canonical) == 8192
    assert all(canonical_rc(kmer) == canonical_rc(reverse_complement(kmer)) for kmer in kmers)
