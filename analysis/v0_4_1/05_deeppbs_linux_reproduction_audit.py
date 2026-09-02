from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
EXT = ensure_dir(ROOT / "external" / "deeppbs")
DEEPPBS = EXT / "DeepPBS"
METADATA = ensure_dir(ROOT / "metadata" / "v0_4_1")
RESULTS = ensure_dir(ROOT / "results" / "v0_4_1")
TABLES = ensure_dir(RESULTS / "tables")
FIGURES = ensure_dir(RESULTS / "figures")
RUNS = ensure_dir(RESULTS / "external_runs" / "deeppbs" / "official_example")
DOCS = ensure_dir(ROOT / "docs" / "v0_4_1")
TODAY = date.today().isoformat()


def run_command(cmd: list[str]) -> tuple[int, str]:
    try:
        proc = subprocess.run(
            cmd,
            cwd=ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=30,
        )
        return proc.returncode, (proc.stdout or "") + (proc.stderr or "")
    except FileNotFoundError as exc:
        return 127, str(exc)
    except subprocess.TimeoutExpired as exc:
        return 124, (exc.stdout or "") + (exc.stderr or "") + "\nTIMEOUT"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(args: list[str]) -> str:
    code, out = run_command(["git", "-C", str(DEEPPBS), *args])
    return out.strip() if code == 0 else "unknown"


def write_container_files(commit: str) -> None:
    (EXT / "Dockerfile").write_text(
        f"""# Linux wrapper for official DeepPBS reproduction.
# Prefer the official Docker image documented by DeepPBS. This wrapper records
# the exact GitHub checkout used in this repository for provenance; it does not
# modify third-party source code.
FROM aricohen/deeppbs:latest

LABEL org.opencontainers.image.source="https://github.com/timkartar/DeepPBS"
LABEL org.opencontainers.image.revision="{commit}"
LABEL org.opencontainers.image.description="DeepPBS official container wrapper for dbp-offtarget-project v0.4.1"

WORKDIR /app
""",
        encoding="utf-8",
    )
    (EXT / "run_official_example_linux.sh").write_text(
        """#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/results/v0_4_1/external_runs/deeppbs/official_example"
mkdir -p "${OUT_DIR}"

# Requires Docker on Linux and the official DeepPBS image:
#   docker pull aricohen/deeppbs:latest
# The input should be replaced by a valid protein-DNA complex PDB/CIF if the
# local DeepPBS checkout does not include an example structure.
docker run --rm \
  -v "${REPO_ROOT}/external/deeppbs/DeepPBS/run/process/pdb:/app/input:ro" \
  -v "${OUT_DIR}:/output" \
  aricohen/deeppbs:latest /app/input/5x6g.pdb -m
""",
        encoding="utf-8",
    )


def main() -> None:
    commit = git_value(["rev-parse", "HEAD"])
    remote = git_value(["remote", "get-url", "origin"])
    docker_code, docker_out = run_command(["docker", "--version"])
    wsl_code, wsl_out = run_command(["wsl", "--status"])
    env_text = (
        f"date: {TODAY}\n"
        f"host_os: Windows PowerShell environment\n"
        f"docker_exit_code: {docker_code}\n{docker_out}\n"
        f"wsl_exit_code: {wsl_code}\n{wsl_out}\n"
    )
    (RUNS / "environment_check.txt").write_text(env_text, encoding="utf-8")
    write_container_files(commit)

    weights = []
    for path in sorted((DEEPPBS / "run" / "output").glob("*/Model.best.tar")):
        weights.append(
            {
                "weight_file": str(path.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
                "source": "bundled in checked-out official DeepPBS repository",
            }
        )
    pd.DataFrame(weights).to_csv(METADATA / "deeppbs_weight_manifest.csv", index=False)

    prov = f"""# DeepPBS Provenance for v0.4.1

Audit date: {TODAY}

- Official repository: https://github.com/timkartar/DeepPBS
- Local path: `external/deeppbs/DeepPBS`
- Local checkout commit: `{commit}`
- Remote: `{remote}`
- Paper DOI from README: `10.1038/s41592-024-02372-w`
- License: BSD 3-Clause, see `external/deeppbs/DeepPBS/LICENSE.txt`
- Official Docker image documented by README: `aricohen/deeppbs:latest`
- Bundled model checkpoints recorded in `metadata/v0_4_1/deeppbs_weight_manifest.csv`

## v0.4.1 Docker/WSL Execution Status

The current host does not provide Docker (`docker --version` exit code {docker_code}) or an installed WSL distribution (`wsl --status` exit code {wsl_code}). Therefore the official Linux preprocessing/example run was not executed on this host in v0.4.1.

This is an environment limitation, not a model result. No DeepPBS designed-DBP performance is reported from this run.
"""
    (EXT / "PROVENANCE.md").write_text(prov, encoding="utf-8")

    (DOCS / "DEEPPBS_TASK_MAPPING.md").write_text(
        """# v0.4.1 DeepPBS Task Mapping

DeepPBS predicts DNA specificity profiles/PWMs from a protein-DNA complex structure after Linux-oriented structural preprocessing. The designed-DBP benchmark target is a per-protein ranking of 7-mer reverse-complement classes by processed uPBM E-score.

Pre-registered mapping for a future runnable Linux execution:

1. Use only author-released or experimental protein-DNA structures listed in `metadata/v0_4_1/deeppbs_structure_manifest.csv`.
2. Run the official DeepPBS preprocessing and ensemble inference without modifying third-party model code.
3. Convert the predicted PWM/profile to candidate 7-mer scores by summing log probabilities across the aligned structure-defined DNA window.
4. Evaluate only per protein with Spearman, NDCG@1%, NDCG@5%, pairwise ranking accuracy, and top-1% experimental recovery.

The PWM-derived score must not be described as PBM affinity, Kd, or calibrated binding probability. If fewer than four designed DBPs have structures and runnable preprocessing, DeepPBS remains a partial diagnostic baseline.
""",
        encoding="utf-8",
    )

    overlap = pd.read_csv(ROOT / "metadata" / "v0_4" / "baseline_data_overlap_audit.csv")
    overlap = overlap[overlap["baseline"] == "DeepPBS"].copy()
    overlap.to_csv(METADATA / "deeppbs_overlap_audit.csv", index=False)
    struct = pd.read_csv(ROOT / "metadata" / "v0_4" / "designed_dbp_structure_manifest.csv")
    struct.to_csv(METADATA / "deeppbs_structure_manifest.csv", index=False)

    pred_cols = ["protein_id", "canonical_7mer", "deeppbs_score", "structure_id", "structure_type", "model_version", "overlap_status"]
    pd.DataFrame(columns=pred_cols).to_parquet(TABLES / "deeppbs_designed_predictions.parquet", index=False)
    perf_rows = []
    for protein_id in ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]:
        perf_rows.append(
            {
                "protein_id": protein_id,
                "baseline": "DeepPBS",
                "spearman": pd.NA,
                "ndcg_1pct": pd.NA,
                "ndcg_5pct": pd.NA,
                "pairwise_accuracy": pd.NA,
                "top1pct_recovery": pd.NA,
                "n_rc_classes": 0,
                "evaluation_status": "not_run_host_lacks_docker_or_wsl",
            }
        )
    pd.DataFrame(perf_rows).to_csv(TABLES / "deeppbs_performance.csv", index=False)

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.bar(["evaluable/run", "not run"], [0, 7], color=["#3A6EA5", "#C95D63"])
    ax.set_ylabel("Designed DBPs")
    ax.set_title("DeepPBS v0.4.1 Execution Coverage")
    ax.text(1, 7, "Docker/WSL unavailable", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_v0_4_1_5_deeppbs_designed_performance.png", dpi=300)
    plt.close(fig)
    print(env_text.encode("ascii", "replace").decode("ascii"))
    print(prov.encode("ascii", "replace").decode("ascii"))


if __name__ == "__main__":
    main()
