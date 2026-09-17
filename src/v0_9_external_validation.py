"""One-shot v0.9 external validation runner with a hard label-access firewall."""
from __future__ import annotations

import argparse
import hashlib
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata/v0_9_external/FROZEN_VALIDATION_MANIFEST.yaml"
CONFIG = ROOT / "metadata/v0_9_external/external_validation_config.yaml"
OUT = ROOT / "results/v0_9_external/PRIMARY_CONFIRMATORY_RESULTS.csv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def current_commit() -> str:
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


def assert_frozen_inputs(manifest: dict) -> None:
    if manifest.get("label_access") is not True:
        raise RuntimeError("External validation locked: metadata/v0_9_external/FROZEN_VALIDATION_MANIFEST.yaml label_access must be true after preregistration and cohort freeze.")
    if manifest.get("git_commit") != current_commit():
        raise RuntimeError("External validation locked: current git commit does not match frozen manifest.")
    for rel, expected in manifest.get("code_hashes", {}).items():
        path = ROOT / rel
        if expected.startswith("PENDING") or not path.exists() or sha256(path) != expected:
            raise RuntimeError(f"External validation locked: code/config hash mismatch for {rel}.")
    if not CONFIG.exists():
        raise RuntimeError("External validation locked: frozen config is missing.")
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    if config.get("label_access") is not True:
        raise RuntimeError("External validation locked: config label_access must also be true.")
    if config.get("status") != "frozen_before_external_label_access":
        raise RuntimeError("External validation locked: unexpected config status.")


def validate_dataset_schema(dataset: Path) -> pd.DataFrame:
    if not dataset.exists():
        raise RuntimeError(f"Dataset is missing: {dataset}")
    frame = pd.read_csv(dataset)
    required = {"study", "assay_tier", "exposure_status", "label_access"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"Dataset schema missing columns: {sorted(missing)}")
    label_flags = frame["label_access"].astype(str).str.strip().str.lower().eq("true")
    if label_flags.any():
        raise RuntimeError("Dataset schema contains label-accessible rows in the pre-unlock manifest.")
    return frame


def validate_label_schema(labels: Path) -> pd.DataFrame:
    if not labels.exists():
        raise RuntimeError(f"External labels are missing: {labels}")
    frame = pd.read_csv(labels)
    required = {"protein_id", "candidate_dna", "measured_score"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"External label schema missing columns: {sorted(missing)}")
    if frame.empty:
        raise RuntimeError("External label table is empty.")
    return frame


def run_once(dataset: Path, labels: Path | None) -> None:
    manifest = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))
    assert_frozen_inputs(manifest)
    audit = validate_dataset_schema(dataset)
    eligible = audit[audit.exposure_status.str.startswith("UNEXPOSED_CONFIRMATORY")].copy()
    if eligible.empty:
        raise RuntimeError("No qualifying external cohort in the frozen dataset manifest.")
    if labels is None:
        raise RuntimeError("External validation locked: pass --labels only after the cohort, schema, and manifest are frozen.")
    validate_label_schema(labels)
    raise RuntimeError("External validation fail-closed: no cohort-specific method adapters are registered; no primary result was written.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=ROOT / "metadata/v0_9_external/external_dataset_audit.csv")
    parser.add_argument("--labels", type=Path, default=None, help="Unlocked external label table with protein_id,candidate_dna,measured_score")
    args = parser.parse_args()
    run_once(args.dataset, args.labels)


if __name__ == "__main__":
    main()
