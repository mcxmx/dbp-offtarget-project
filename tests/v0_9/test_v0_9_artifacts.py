from pathlib import Path
import subprocess

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def test_v09_required_files_exist():
    required = [
        "metadata/v0_9_external/external_dataset_audit.csv",
        "metadata/v0_9_external/external_validation_config.yaml",
        "metadata/v0_9_external/FROZEN_VALIDATION_MANIFEST.yaml",
        "docs/v0_9_external/ASSAY_SPECIFIC_VALIDATION_RULES.md",
        "docs/v0_9_external/RFD3_COHORT_AUDIT.md",
        "docs/v0_9_external/DATA_REQUEST_EMAIL_RFD3.md",
        "docs/v0_9_external/PREREGISTERED_EXTERNAL_VALIDATION.md",
        "docs/paper/MANUSCRIPT_V0_9.md",
        "docs/paper/REVIEWER_ATTACKS.md",
        "docs/paper/SUBMISSION_READINESS.md",
        "docs/paper/FIGURE_PLAN_V0_9.md",
        "src/v0_9_external_validation.py",
    ]
    assert all((ROOT / path).exists() for path in required)


def test_external_audit_is_label_locked():
    frame = pd.read_csv(ROOT / "metadata/v0_9_external/external_dataset_audit.csv")
    assert frame["label_access"].astype(str).str.lower().eq("false").all()
    rfd3 = frame.loc[frame.study.str.startswith("RFdiffusion3")].iloc[0]
    assert rfd3.assay_tier == "Tier 2 candidate"
    assert rfd3.exposure_status == "UNEXPOSED_CONFIRMATORY_CANDIDATE"
    manifest = (ROOT / "metadata/v0_9_external/FROZEN_VALIDATION_MANIFEST.yaml").read_text(encoding="utf-8")
    assert "PENDING" not in manifest


def test_runner_refuses_locked_state_and_does_not_write_results(tmp_path):
    output = ROOT / "results/v0_9_external/PRIMARY_CONFIRMATORY_RESULTS.csv"
    if output.exists():
        output.unlink()
    proc = subprocess.run(
        [str(ROOT / ".venv313/Scripts/python.exe"), "src/v0_9_external_validation.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    assert "locked" in (proc.stdout + proc.stderr).lower()
    assert not output.exists()


def test_manuscript_claim_discipline():
    manuscript = (ROOT / "docs/paper/MANUSCRIPT_V0_9.md").read_text(encoding="utf-8")
    assert "sequence-based" in manuscript
    assert "universal claim" in manuscript
    assert "Reserved for preregistered independent validation" in manuscript
    forbidden = ["computational methods fail", "structure-based methods fail", "experimental upper bound", "genome-wide off-target risk predictor"]
    assert not any(term in manuscript for term in forbidden)
