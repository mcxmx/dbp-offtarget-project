from hashlib import sha256
from pathlib import Path

import pandas as pd

from src.v0_6.gates import summarize_replay_gates


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "v0_6_representation"


def test_v06_required_machine_readable_artifacts_exist_and_have_expected_rows():
    required = {
        "representation_audit.csv": 1,
        "synthetic_conditionality.csv": 4,
        "primary_results.csv": 112,
        "protein_shuffle.csv": 21,
        "target_shuffle.csv": 21,
        "conditionality_diagnostics.csv": 5,
        "go_no_go.csv": 10,
    }
    for filename, rows in required.items():
        frame = pd.read_csv(RESULTS / filename)
        assert len(frame) == rows, filename


def test_v06_required_reports_state_developmental_exposure():
    audit = (ROOT / "docs" / "v0_6_representation" / "V0_6_REPRESENTATION_AUDIT.md").read_text(encoding="utf-8")
    results = (ROOT / "docs" / "v0_6_representation" / "V0_6_RESULTS.md").read_text(encoding="utf-8")
    assert "not independent confirmatory" in audit
    assert "seven designed DBPs are all development-exposed" in results
    exposure = pd.read_csv(ROOT / "metadata" / "v0_6" / "exposure_manifest.csv")
    assert exposure["independent_confirmatory"].eq(False).all()
    assert exposure["v0_6_use"].eq("developmental_debugging_only").all()


def test_v05_frozen_result_hashes_remain_unchanged():
    manifest = ROOT / "results" / "v0_5" / "PRIMARY_RESULTS_FROZEN_MANIFEST.txt"
    for line in manifest.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) != 2 or not parts[0].startswith("results/v0_5/"):
            continue
        path, expected = parts
        assert sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
