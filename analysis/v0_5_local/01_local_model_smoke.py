from __future__ import annotations

from src.v0_5_local_training import LocalConfig, run_local_smoke, write_local_smoke_outputs


def main() -> None:
    config = LocalConfig.from_json()
    outputs = run_local_smoke(config)
    write_local_smoke_outputs(outputs)
    print(outputs["smoke_info"].to_string(index=False))
    print(outputs["smoke_results"].to_string(index=False))
    print(outputs["health"].to_string(index=False))


if __name__ == "__main__":
    main()
