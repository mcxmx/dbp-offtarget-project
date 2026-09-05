from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root
from src.v0_5_primary_evaluation import build_primary_artifacts
from src.v0_5_training import V05Config


def write_artifacts(artifacts: dict[str, pd.DataFrame], config: V05Config) -> None:
    root = project_root()
    results = ensure_dir(root / "results" / "v0_5")
    logs = ensure_dir(root / "logs" / "v0_5")
    metadata = ensure_dir(root / "metadata" / "v0_5")
    artifacts["primary_seed"].to_csv(results / "primary_seed_level_results.csv", index=False)
    artifacts["primary_per_protein"].to_csv(results / "primary_per_protein_results.csv", index=False)
    artifacts["primary_macro"].to_csv(results / "primary_macro_summary.csv", index=False)
    artifacts["strict_seed"].to_csv(results / "strict_component_seed_level_results.csv", index=False)
    artifacts["strict_per_protein"].to_csv(results / "strict_component_per_protein_results.csv", index=False)
    artifacts["strict_macro"].to_csv(results / "strict_component_macro_summary.csv", index=False)
    artifacts["controls"].to_csv(results / "target_relative_controls_full.csv", index=False)
    artifacts["baseline_context"].to_csv(results / "baseline_context_table.csv", index=False)
    artifacts["training_health"].to_csv(results / "primary_training_health.csv", index=False)
    (metadata / "v0_5_model_config.json").write_text(
        json.dumps(config.as_dict(), indent=2),
        encoding="utf-8",
    )
    log_lines = [
        "v0.5 complete primary and strict evaluation",
        f"evaluation_seeds={','.join(str(seed) for seed in config.evaluation_seeds)}",
        artifacts["primary_per_protein"].to_string(index=False),
        artifacts["strict_macro"].to_string(index=False),
        artifacts["controls"].to_string(index=False),
    ]
    (logs / "primary_evaluation.log").write_text("\n\n".join(log_lines), encoding="utf-8")


def main() -> None:
    config = V05Config()
    artifacts = build_primary_artifacts(config, project_root())
    write_artifacts(artifacts, config)
    print(artifacts["primary_macro"].to_string(index=False))
    print(artifacts["strict_macro"].to_string(index=False))


if __name__ == "__main__":
    main()
