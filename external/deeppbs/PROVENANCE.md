# DeepPBS Provenance for v0.4.1

Audit date: 2026-09-02

- Official repository: https://github.com/timkartar/DeepPBS
- Local path: `external/deeppbs/DeepPBS`
- Local checkout commit: `8bfb211dd67f02877841f6f33aa493ddf7daedf9`
- Remote: `https://github.com/timkartar/DeepPBS.git`
- Paper DOI from README: `10.1038/s41592-024-02372-w`
- License: BSD 3-Clause, see `external/deeppbs/DeepPBS/LICENSE.txt`
- Official Docker image documented by README: `aricohen/deeppbs:latest`
- Bundled model checkpoints recorded in `metadata/v0_4_1/deeppbs_weight_manifest.csv`

## v0.4.1 Docker/WSL Execution Status

The current host does not provide Docker (`docker --version` exit code 127) or an installed WSL distribution (`wsl --status` exit code 50). Therefore the official Linux preprocessing/example run was not executed on this host in v0.4.1.

This is an environment limitation, not a model result. No DeepPBS designed-DBP performance is reported from this run.
