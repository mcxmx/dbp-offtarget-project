from __future__ import annotations

from .audit import run_representation_audit
from .gates import assert_pre_replay_gates
from .report import generate_reports
from .replay import run_v06_replay


def main() -> None:
    audit = run_representation_audit()
    gates = assert_pre_replay_gates()
    gates.to_csv("results/v0_6_representation/pre_replay_gates.csv", index=False)
    print("v0.6 pre-replay gates passed")
    for frame in audit.values():
        print(frame.to_string(index=False))
    outputs = run_v06_replay()
    outputs["health"].to_csv("results/v0_6_representation/training_health.csv", index=False)
    generate_reports()
    print(outputs["primary"].head().to_string(index=False))


if __name__ == "__main__":
    main()
