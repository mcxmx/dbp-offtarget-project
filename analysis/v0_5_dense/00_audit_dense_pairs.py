from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.utils import ensure_dir, project_root
from src.v0_5_dense_sampling import audit_pair_sampling
from src.v0_5_training import load_v05_data


ROOT = project_root()
OUTPUT_DIR = ROOT / "results" / "v0_5_dense"


def main() -> None:
    ensure_dir(OUTPUT_DIR)
    benchmark, _, _, _ = load_v05_data()
    proteins = sorted(benchmark["protein_id"].unique())
    audit = audit_pair_sampling(benchmark, proteins, seed=17)
    audit.to_csv(OUTPUT_DIR / "pair_sampling_audit.csv", index=False)
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
