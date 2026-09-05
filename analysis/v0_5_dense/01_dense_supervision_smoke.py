from __future__ import annotations

from src.v0_5_dense_training import run_dense_smoke, write_dense_outputs


def main() -> None:
    outputs = run_dense_smoke()
    write_dense_outputs(outputs)
    print(outputs["info"].to_string(index=False))
    print(outputs["health"].to_string(index=False))
    print(outputs["delta"].to_string(index=False))
    print(outputs["shuffle"].to_string(index=False))


if __name__ == "__main__":
    main()
