from __future__ import annotations

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.v0_5_training import V05Config, run_smoke, write_smoke_outputs


def main() -> None:
    config = V05Config()
    outputs = run_smoke(config)
    write_smoke_outputs(config, outputs)
    print(outputs["smoke_info"].to_string(index=False))
    print(outputs["evaluation"].to_string(index=False))
    print(outputs["controls"].to_string(index=False))
    print(outputs["parameters"].to_string(index=False))


if __name__ == "__main__":
    main()
