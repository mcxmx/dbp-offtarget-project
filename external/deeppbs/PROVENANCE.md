# DeepPBS Provenance for v0.4.2

Audit date: 2026-09-04

- Official repository: https://github.com/timkartar/DeepPBS
- Local checkout: `external/deeppbs/DeepPBS`
- Local checkout commit: `8bfb211dd67f02877841f6f33aa493ddf7daedf9`
- Remote: `https://github.com/timkartar/DeepPBS.git`
- Paper DOI from official README: `10.1038/s41592-024-02372-w`
- License: BSD 3-Clause, see `external/deeppbs/DeepPBS/LICENSE.txt`
- Official Docker image documented by README: `aricohen/deeppbs:latest`
- Bundled model checkpoints: `metadata/v0_4_2/deeppbs_weight_manifest_v2.csv`
- Completed execution host: Ubuntu 20.04 VM, `qwqaq@192.168.73.128`
- Python: `/home/qwqaq/miniconda3/envs/deeppbs/bin/python`
- Python version: 3.10.14
- PyTorch: 2.2.1+cpu
- NumPy: 1.26.3
- SciPy: 1.12.0
- Biopython: 1.83
- X3DNA: upstream `x3dna-v2.3-linux-64bit`
- GPU: none; CPU inference

## Execution status

The official `5x6g` example passed in the Ubuntu VM. Its logs and NPZ
artifacts are retained under
`results/v0_4_2/external_runs/deeppbs_official_example/`, with the acceptance
marker `OFFICIAL_EXAMPLE_PASS`.

The project then ran all legally evaluable designed DBP inputs:

- DBP35: original project design complex, return code 0, Helix score 1.0,
  contact count 266.
- DBP48: PDB 8TAC after a DSSR-defined 9-bp helix-only input preparation,
  return code 0, Helix score 1.0, contact count 224.

The original 8TAC file is preserved. The upstream model and weights were not
modified. Earlier Docker/WSL execution was host-limited; the completed run
used the Ubuntu VM instead.

Exact commands, stdout/stderr, input/output paths, and checksums are recorded
in `results/v0_4_2/tables/deeppbs_run_manifest_completed_v0_4_2.csv`.
The source-level interpretation of `P` and `Seq` is in
`docs/v0_4_2/DEEPPBS_OUTPUT_SEMANTICS.md`.
