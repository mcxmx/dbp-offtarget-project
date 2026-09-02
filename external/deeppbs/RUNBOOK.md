# DeepPBS v0.4.2 Linux Runbook

The official DeepPBS preprocessing is Linux-oriented and requires a Linux container/runtime. This Windows host has no Docker, Podman, conda/mamba, or installed WSL distribution, so the official example was not executed locally.

## Primary official Docker path

```bash
docker pull aricohen/deeppbs:latest
docker run --rm \
  -v "$PWD/external/deeppbs/DeepPBS/run/process/pdb:/app/input:ro" \
  -v "$PWD/results/v0_4_2/external_runs/deeppbs_official_example:/output" \
  aricohen/deeppbs:latest /app/input/5x6g.pdb -m
```

Expected first acceptance artifact: `results/v0_4_2/external_runs/deeppbs_official_example/predict/`.

## Local checkout

- Repository: https://github.com/timkartar/DeepPBS.git
- Commit: 8bfb211dd67f02877841f6f33aa493ddf7daedf9
- License: BSD 3-Clause
- Weights: `metadata/v0_4_2/deeppbs_weight_manifest_v2.csv`

No third-party DeepPBS source code is modified by this project.
