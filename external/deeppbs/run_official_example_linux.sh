#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/results/v0_4_1/external_runs/deeppbs/official_example"
mkdir -p "${OUT_DIR}"

# Requires Docker on Linux and the official DeepPBS image:
#   docker pull aricohen/deeppbs:latest
# The input should be replaced by a valid protein-DNA complex PDB/CIF if the
# local DeepPBS checkout does not include an example structure.
docker run --rm   -v "${REPO_ROOT}/external/deeppbs/DeepPBS/run/process/pdb:/app/input:ro"   -v "${OUT_DIR}:/output"   aricohen/deeppbs:latest /app/input/5x6g.pdb -m
