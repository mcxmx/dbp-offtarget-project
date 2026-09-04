# DeepPBS v0.4.2 Linux Runbook

The official DeepPBS preprocessing is Linux-oriented. The reproducible
execution used the configured Ubuntu VM rather than native Windows.

## Completed environment

- SSH host: `qwqaq@192.168.73.128`
- Upstream checkout: `/home/qwqaq/DeepPBS`
- Environment Python: `/home/qwqaq/miniconda3/envs/deeppbs/bin/python`
- Upstream commit: `8bfb211dd67f02877841f6f33aa493ddf7daedf9`
- Execution mode: CPU

The official example was run first and passed. The captured acceptance marker
is `results/v0_4_2/external_runs/deeppbs_official_example/official_example_status.txt`.

## Documented upstream command pattern

The following is the command pattern used remotely, with the actual per-input
filenames recorded in each project-side run manifest:

```bash
export PATH=/home/qwqaq/miniconda3/envs/deeppbs/bin:/home/qwqaq/DeepPBS/dependencies/bin:$PATH
export X3DNA=/home/qwqaq/DeepPBS/x3dna-v2.3-linux-64bit/x3dna-v2.3
cd /home/qwqaq/DeepPBS/run/process
/home/qwqaq/miniconda3/envs/deeppbs/bin/python ../process_co_crystal.py INPUT.pdb CONFIG --no_pwm
/home/qwqaq/miniconda3/envs/deeppbs/bin/python ../predict.py INPUT.npz OUTPUT_DIR -c CONFIG
```

DBP48 required a project-side DSSR-defined helix-only PDB because the original
deposited DNA contains non-helical overhang residues that fail upstream shape
extraction. This preparation is documented and does not modify the upstream
checkout.

## Provenance

- Repository: https://github.com/timkartar/DeepPBS.git
- License: BSD 3-Clause
- Weights: `metadata/v0_4_2/deeppbs_weight_manifest_v2.csv`
- Full run manifest: `results/v0_4_2/tables/deeppbs_run_manifest_completed_v0_4_2.csv`
- Output semantics: `docs/v0_4_2/DEEPPBS_OUTPUT_SEMANTICS.md`

No third-party DeepPBS source code is modified by this project.
