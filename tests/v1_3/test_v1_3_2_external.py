from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.v1_3.run_v1_3_2_external import assign_split, pairwise_accuracy


ROOT = Path(__file__).resolve().parents[2]


def test_primary_external_outputs_have_frozen_strata_and_counts() -> None:
    spots = pd.read_parquet(ROOT / "data/processed/v1_3_2_samba_spots.parquet")
    assert len(spots) == 22_026
    assert spots.tf_id.nunique() == 7
    assert spots.binding_site_id.nunique() == 12
    assert set(spots.perturbation_type) == {
        "reference",
        "canonical_watson_crick",
        "mismatch",
    }
    assert set(spots.array_id) == {"unknown_array"}
    assert set(spots.array_id_source) == {"not_supplied_by_source_workbook"}

    decomposition = pd.read_csv(
        ROOT / "results/v1_3/external_samba_decomposition.tsv", sep="\t"
    )
    assert set(decomposition.perturbation_type) == {
        "canonical_watson_crick",
        "mismatch",
    }
    assert len(decomposition) == 20
    assert (decomposition.n_positions >= 3).all()


def test_technical_split_is_deterministic_and_balanced() -> None:
    frame = pd.DataFrame(
        {
            "tf_id": ["T"] * 6,
            "binding_site_id": ["S"] * 6,
            "sequence": ["A"] * 3 + ["B"] * 3,
            "array_id": ["unknown_array"] * 6,
            "spot_id": [f"spot_{i}" for i in range(6)],
            "raw_signal": np.arange(1, 7, dtype=float),
            "source_mapping_status": ["VALIDATED_BY_EXACT_MEDIAN_RECONSTRUCTION"] * 6,
        }
    )
    first = assign_split(frame, 1301)
    second = assign_split(frame, 1301)
    pd.testing.assert_series_equal(first.technical_split, second.technical_split)
    counts = first.groupby("sequence").technical_split.value_counts()
    assert set(counts.index.get_level_values("technical_split")) == {"A", "B"}
    for _, sequence_counts in counts.groupby(level="sequence"):
        assert sequence_counts.max() - sequence_counts.min() <= 1


def test_pairwise_tie_handling_matches_contract() -> None:
    score, n_pairs = pairwise_accuracy(np.array([2.0, 1.0, 0.0]), np.array([1.0, 1.0, 0.0]))
    assert n_pairs == 3
    assert score == (0.5 + 1.0 + 1.0) / 3


def test_external_manifest_contains_all_downloaded_samba_sources() -> None:
    manifest = json.loads(
        (ROOT / "artifacts/cache/file_manifest.json").read_text(encoding="utf-8")
    )
    paths = set(manifest["files"])
    expected = {
        "data/raw/v1_3_2_external/samba/41586_2020_2843_MOESM4_ESM.xlsx",
        "data/raw/v1_3_2_external/samba/41586_2020_2843_MOESM7_ESM.xlsx",
        "data/raw/v1_3_2_external/samba/GSE156375_RAW.tar",
    }
    assert expected <= paths
