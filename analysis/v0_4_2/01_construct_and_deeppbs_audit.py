from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
TODAY = date.today().isoformat()
V041_META = ROOT / "metadata" / "v0_4_1"
V04_META = ROOT / "metadata" / "v0_4"
V042_META = ensure_dir(ROOT / "metadata" / "v0_4_2")
V042_DOCS = ensure_dir(ROOT / "docs" / "v0_4_2")
V042_TABLES = ensure_dir(ROOT / "results" / "v0_4_2" / "tables")
V042_RUNS = ensure_dir(ROOT / "results" / "v0_4_2" / "external_runs")
V042_PROCESSED = ensure_dir(ROOT / "data" / "processed" / "v0_4_2")
DEEPPBS_EXT = ensure_dir(ROOT / "external" / "deeppbs")
DEEPPBS_REPO = DEEPPBS_EXT / "DeepPBS"


def run_command(cmd: list[str], cwd: Path | None = None, timeout: int = 60) -> dict[str, str | int]:
    try:
        proc = subprocess.run(
            cmd,
            cwd=cwd or ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
        return {"returncode": proc.returncode, "stdout": proc.stdout or "", "stderr": proc.stderr or ""}
    except FileNotFoundError as exc:
        return {"returncode": 127, "stdout": "", "stderr": str(exc)}
    except subprocess.TimeoutExpired as exc:
        return {
            "returncode": 124,
            "stdout": exc.stdout or "",
            "stderr": (exc.stderr or "") + "\nTIMEOUT",
        }


def git_value(args: list[str], repo: Path = DEEPPBS_REPO) -> str:
    result = run_command(["git", "-C", str(repo), *args])
    if result["returncode"] == 0:
        return str(result["stdout"]).strip()
    return "unknown"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_construct_audit() -> None:
    natural_benchmark = pd.read_parquet(
        ROOT / "data" / "processed" / "v0_4_1" / "natural_pbm_benchmark_v0_4_1.parquet"
    )
    benchmark_ids = set(natural_benchmark["protein_id"].drop_duplicates())
    proteins = pd.read_csv(V041_META / "natural_pbm_proteins.csv")
    proteins = proteins[proteins["natural_protein_id"].isin(benchmark_ids)].copy()
    unique = (
        proteins.sort_values(["natural_protein_id", "experiment_id"])
        .drop_duplicates("natural_protein_id")
        .rename(columns={"natural_protein_id": "protein_id"})
        .copy()
    )
    audit_rows = []
    for _, row in unique.iterrows():
        seq_match = str(row.get("sequence_match_to_assay", "")).lower() == "true"
        if seq_match:
            construct_type = row.get("sequence_type", "unknown")
            construct_sequence = row.get("protein_sequence", pd.NA)
            construct_known = True
            confidence = "medium"
            use_primary = True
            notes = "v0.4.1 metadata marked sequence_match_to_assay=True; no additional manual override applied."
        else:
            construct_type = "unknown"
            construct_sequence = pd.NA
            construct_known = False
            confidence = "low_construct_unknown"
            use_primary = False
            notes = (
                "UniPROBE v0.4.1 benchmark uses a UniProt/reference full-length sequence. "
                "The experimental PBM construct sequence, domain boundaries, fusion tag, "
                "and isoform are not recovered in the current local provenance."
            )
        audit_rows.append(
            {
                "protein_id": row["protein_id"],
                "experiment_id": row.get("experiment_id", ""),
                "reference_full_length_sequence": row.get("protein_sequence", ""),
                "experimental_construct_known": construct_known,
                "construct_sequence": construct_sequence,
                "construct_start": pd.NA,
                "construct_end": pd.NA,
                "construct_type": construct_type,
                "construct_source": "v0.4.1 UniPROBE/UniProt metadata audit",
                "fusion_or_tag": "unknown",
                "sequence_confidence": confidence,
                "use_for_primary_training": use_primary,
                "notes": notes,
            }
        )
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(V042_META / "natural_pbm_construct_audit.csv", index=False)

    exact = audit["construct_type"].eq("experimental_construct").sum()
    reconstructed = audit["construct_type"].isin(["dna_binding_domain", "truncated"]).sum()
    confirmed_full = (audit["construct_type"].eq("full_length") & audit["experimental_construct_known"]).sum()
    unknown = audit["construct_type"].eq("unknown").sum()
    high_conf = audit["use_for_primary_training"].sum()
    coverage = pd.DataFrame(
        [
            {
                "n_total_proteins": int(len(audit)),
                "n_exact_construct": int(exact),
                "n_domain_reconstructed": int(reconstructed),
                "n_full_length_confirmed": int(confirmed_full),
                "n_unknown": int(unknown),
                "fraction_high_confidence_construct": float(high_conf / len(audit)) if len(audit) else 0.0,
                "interpretation": (
                    "Current v0.4.2 audit cannot verify assay construct sequences for the "
                    "natural UniPROBE proteins; assay-aligned primary training is therefore empty."
                ),
            }
        ]
    )
    coverage.to_csv(V042_TABLES / "natural_construct_coverage.csv", index=False)

    full_length = natural_benchmark.copy()
    full_length["protein_sequence_version"] = "FULL_LENGTH_REFERENCE"
    full_length.to_parquet(V042_PROCESSED / "natural_pbm_full_length_reference_v0_4_2.parquet", index=False)

    assay_aligned_cols = list(natural_benchmark.columns) + [
        "assay_aligned_sequence",
        "assay_aligned_sequence_type",
        "assay_aligned_sequence_source",
    ]
    assay_aligned = pd.DataFrame(columns=assay_aligned_cols)
    assay_aligned.to_parquet(V042_PROCESSED / "natural_pbm_assay_aligned_v0_4_2.parquet", index=False)

    splits = pd.read_csv(V041_META / "natural_pbm_splits.csv").rename(columns={"protein_id": "protein_id"})
    splits.to_csv(V042_META / "natural_pbm_full_length_reference_splits.csv", index=False)
    pd.DataFrame(columns=list(splits.columns) + ["assay_aligned_status"]).to_csv(
        V042_META / "natural_pbm_construct_aware_splits.csv", index=False
    )

    (V042_DOCS / "NATURAL_PBM_CONSTRUCT_AUDIT.md").write_text(
        f"""# v0.4.2 Natural PBM Construct Audit

Audit date: {TODAY}

The v0.4.1 natural PBM benchmark recovered protein sequences primarily from UniProt/reference records. The local metadata explicitly does not claim those sequences are the exact PBM assay constructs.

## Result

- Natural proteins audited: {len(audit)}
- Exact experimental constructs recovered: {int(exact)}
- Domain/truncated constructs reconstructed from reported coordinates: {int(reconstructed)}
- Confirmed full-length assay constructs: {int(confirmed_full)}
- Unknown assay constructs: {int(unknown)}
- High-confidence construct coverage: {float(high_conf / len(audit)) if len(audit) else 0.0:.3f}

## Consequence

`FULL_LENGTH_REFERENCE` remains a sensitivity benchmark because it is reproducible and has sequence provenance, but it is not an assay-aligned construct benchmark. `ASSAY_ALIGNED_PROTEIN` is empty in v0.4.2 because no construct sequence or reliable construct coordinate provenance was recovered.

No missing construct sequence is filled by guessing, domain heuristics, or family-level substitution.
""",
        encoding="utf-8",
    )


def write_deeppbs_audit() -> None:
    commit = git_value(["rev-parse", "HEAD"])
    remote = git_value(["remote", "get-url", "origin"])
    license_text = "BSD 3-Clause" if (DEEPPBS_REPO / "LICENSE.txt").exists() else "unknown"
    weights = []
    for path in sorted((DEEPPBS_REPO / "run" / "output").glob("*/Model.best.tar")):
        weights.append(
            {
                "weight_file": str(path.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
                "source": "bundled in checked-out official DeepPBS repository",
            }
        )
    pd.DataFrame(weights).to_csv(V042_META / "deeppbs_weight_manifest_v2.csv", index=False)

    docker = run_command(["docker", "--version"])
    podman = run_command(["podman", "--version"])
    wsl = run_command(["wsl", "--status"])
    runtime_available = docker["returncode"] == 0 or podman["returncode"] == 0 or wsl["returncode"] == 0

    example_dir = ensure_dir(V042_RUNS / "deeppbs_official_example")
    command = (
        "docker run --rm -v ${PWD}/external/deeppbs/DeepPBS/run/process/pdb:/app/input:ro "
        "-v ${PWD}/results/v0_4_2/external_runs/deeppbs_official_example:/output "
        "aricohen/deeppbs:latest /app/input/5x6g.pdb -m"
    )
    (example_dir / "command.txt").write_text(command + "\n", encoding="utf-8")
    (example_dir / "stdout.txt").write_text("", encoding="utf-8")
    (example_dir / "stderr.txt").write_text(
        "OFFICIAL_EXAMPLE_NOT_RUN_HOST_NO_LINUX_RUNTIME\n"
        f"docker_returncode={docker['returncode']}\n{docker['stderr']}\n"
        f"podman_returncode={podman['returncode']}\n{podman['stderr']}\n"
        f"wsl_returncode={wsl['returncode']}\n{wsl['stdout']}{wsl['stderr']}\n",
        encoding="utf-8",
    )
    (example_dir / "runtime.txt").write_text("not_run_host_no_linux_container_runtime\n", encoding="utf-8")
    (example_dir / "environment.txt").write_text(
        json.dumps(
            {
                "audit_date": TODAY,
                "host": "Windows PowerShell",
                "docker_returncode": docker["returncode"],
                "podman_returncode": podman["returncode"],
                "wsl_returncode": wsl["returncode"],
                "official_example_status": "OFFICIAL_EXAMPLE_NOT_RUN_HOST_NO_LINUX_RUNTIME",
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    (DEEPPBS_EXT / "Dockerfile").write_text(
        f"""# Linux wrapper for official DeepPBS reproduction.
FROM aricohen/deeppbs:latest

LABEL org.opencontainers.image.source="https://github.com/timkartar/DeepPBS"
LABEL org.opencontainers.image.revision="{commit}"
LABEL org.opencontainers.image.description="DeepPBS official container wrapper for dbp-offtarget-project v0.4.2"

WORKDIR /app
""",
        encoding="utf-8",
    )
    (DEEPPBS_EXT / "environment.yml").write_text(
        """name: deeppbs-linux-official
channels:
  - pytorch
  - nvidia
  - pyg
  - conda-forge
dependencies:
  - python=3.10
  - pytorch=2.3.0
  - torchvision=0.18.0
  - torchaudio=2.3.0
  - pytorch-cuda=12.1
  - pip
  - pip:
      - torch_geometric
      - biopython==1.83
      - logomaker
      - matplotlib==3.5.2
      - networkx
      - pandas==1.4.4
      - pdb2pqr
      - scipy==1.14.1
      - seaborn==0.13.2
      - freesasa==2.2.1
""",
        encoding="utf-8",
    )
    (DEEPPBS_EXT / "RUNBOOK.md").write_text(
        f"""# DeepPBS v0.4.2 Linux Runbook

The official DeepPBS preprocessing is Linux-oriented and requires a Linux container/runtime. This Windows host has no Docker, Podman, conda/mamba, or installed WSL distribution, so the official example was not executed locally.

## Primary official Docker path

```bash
docker pull aricohen/deeppbs:latest
docker run --rm \\
  -v "$PWD/external/deeppbs/DeepPBS/run/process/pdb:/app/input:ro" \\
  -v "$PWD/results/v0_4_2/external_runs/deeppbs_official_example:/output" \\
  aricohen/deeppbs:latest /app/input/5x6g.pdb -m
```

Expected first acceptance artifact: `results/v0_4_2/external_runs/deeppbs_official_example/predict/`.

## Local checkout

- Repository: {remote}
- Commit: {commit}
- License: {license_text}
- Weights: `metadata/v0_4_2/deeppbs_weight_manifest_v2.csv`

No third-party DeepPBS source code is modified by this project.
""",
        encoding="utf-8",
    )
    (DEEPPBS_EXT / "PROVENANCE.md").write_text(
        f"""# DeepPBS Provenance for v0.4.2

Audit date: {TODAY}

- Official repository: https://github.com/timkartar/DeepPBS
- Local checkout: `external/deeppbs/DeepPBS`
- Local checkout commit: `{commit}`
- Remote: `{remote}`
- Paper DOI from official README: `10.1038/s41592-024-02372-w`
- License: {license_text}, see `external/deeppbs/DeepPBS/LICENSE.txt`
- Official Docker image documented by README: `aricohen/deeppbs:latest`
- Bundled model checkpoints: `metadata/v0_4_2/deeppbs_weight_manifest_v2.csv`

## v0.4.2 Execution Status

Official example status: `OFFICIAL_EXAMPLE_NOT_RUN_HOST_NO_LINUX_RUNTIME`.

This is an environment/runtime limitation, not a DeepPBS model performance result. No DeepPBS designed-DBP ranking is reported from this host.
""",
        encoding="utf-8",
    )

    (V042_DOCS / "DEEPPBS_SCORING_PROTOCOL.md").write_text(
        """# DeepPBS v0.4.2 Scoring Protocol

This protocol is frozen before any designed-DBP DeepPBS result is evaluated.

## Inputs

Use only public/original protein-DNA complex structures listed in `metadata/v0_4_2/designed_structure_manifest_v2.csv`.

## Primary mapping

If official DeepPBS returns a position-wise PWM/profile over a structure-defined DNA binding window, the primary candidate 7-mer score is:

`score(D) = max_offset sum_i log(P(D_i at aligned position offset+i) + 1e-9)`

The offset range is all contiguous 7-bp windows fully contained in the structure-defined DeepPBS output window. Reverse-complement candidate classes are scored by the maximum over both orientations.

## Constraints

- Do not tune offsets, checkpoint, or aggregation after inspecting PBM E-scores.
- Do not call the PWM-derived value a PBM affinity, Kd, or calibrated binding probability.
- If fewer than 4 designed DBPs are evaluable, DeepPBS remains diagnostic.
""",
        encoding="utf-8",
    )

    if (V04_META / "baseline_data_overlap_audit.csv").exists():
        overlap = pd.read_csv(V04_META / "baseline_data_overlap_audit.csv")
        overlap = overlap[overlap["baseline"].eq("DeepPBS")].copy()
        out = pd.DataFrame(
            {
                "protein_id": overlap["protein_id"],
                "exact_sequence_seen": overlap["exact_protein_seen"],
                "homolog_seen": overlap["homolog_seen"],
                "structure_seen": overlap["structure_seen"],
                "paper_evaluation_seen": False,
                "training_status": "not_confirmed_seen; homolog audit unresolved",
                "evidence": overlap["evidence"],
                "risk_level": overlap["risk_level"],
                "notes": overlap["notes"],
            }
        )
    else:
        out = pd.DataFrame()
    out.to_csv(V042_META / "deeppbs_overlap_audit_v2.csv", index=False)

    if (V04_META / "designed_dbp_structure_manifest.csv").exists():
        struct = pd.read_csv(V04_META / "designed_dbp_structure_manifest.csv")
        struct_v2 = pd.DataFrame(
            {
                "protein_id": struct["protein_id"],
                "structure_available": ~struct["pdb_id"].fillna("").eq("")
                | struct["local_file"].fillna("").str.endswith(".pdb"),
                "structure_source": struct["source_url"],
                "structure_type": struct["structure_type"],
                "pdb_id": struct["pdb_id"],
                "bound_dna": struct["bound_dna_sequence"],
                "bound_dna_length": struct["bound_dna_length"],
                "model_source": struct["local_file"],
                "paper_reference": struct["paper_reference"],
                "quality_notes": struct["notes"],
                "deeppbs_evaluable": False,
            }
        )
        struct_v2.loc[
            struct_v2["protein_id"].isin(["DBP35", "DBP48"]) & struct_v2["structure_available"],
            "deeppbs_evaluable",
        ] = runtime_available
    else:
        struct_v2 = pd.DataFrame()
    struct_v2.to_csv(V042_META / "designed_structure_manifest_v2.csv", index=False)

    pred_cols = [
        "protein_id",
        "canonical_7mer",
        "deeppbs_score",
        "structure_id",
        "structure_source",
        "alignment_rule",
        "overlap_status",
        "model_commit",
        "weight_id",
    ]
    pd.DataFrame(columns=pred_cols).to_parquet(V042_TABLES / "deeppbs_predictions_v0_4_2.parquet", index=False)
    perf = []
    for protein_id in ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]:
        status = "not_evaluable_no_public_structure"
        if protein_id in {"DBP35", "DBP48"}:
            status = "not_run_host_no_linux_container_runtime"
        perf.append(
            {
                "protein_id": protein_id,
                "baseline": "DeepPBS",
                "spearman": pd.NA,
                "ndcg_1pct": pd.NA,
                "ndcg_5pct": pd.NA,
                "pairwise_accuracy": pd.NA,
                "top1pct_recovery": pd.NA,
                "n_rc_classes": 0,
                "evaluation_status": status,
            }
        )
    pd.DataFrame(perf).to_csv(V042_TABLES / "deeppbs_performance_v0_4_2.csv", index=False)


def main() -> None:
    write_construct_audit()
    write_deeppbs_audit()
    print("wrote v0.4.2 construct and DeepPBS audit artifacts")


if __name__ == "__main__":
    main()
