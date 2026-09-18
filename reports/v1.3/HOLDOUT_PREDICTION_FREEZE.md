# Holdout Prediction Freeze

Freeze status: **PREDICTIONS FROZEN; LABELS UNREAD**. Date: 2026-09-17.

This file records Stage 1 of the two-stage prospective evaluation. `scripts/v1_3/freeze_holdout_predictions.py` accepts only DeepPBS NPZ predictions and design PDBs. It does not import, open, or join any competition workbook or processed competition table.

## Frozen artifact

- Prediction file: `results/v1_3/holdout_predictions_unscored.tsv`
- Rows: 129 across 3 proteins
- SHA256: `002c8825679a4a918bf7eac04c06eb390381061a52ae4f5197f5416cbb8dacdf` (computed once and cached in `artifacts/cache/file_manifest.json`)
- Method: official DeepPBS checkout commit `8bfb211dd67f02877841f6f33aa493ddf7daedf9`, bundled five-checkpoint `DeepPBS.txt` ensemble
- Score: `log(P_mutant + 1e-12) - log(P_wt + 1e-12)`; sign unchanged
- Mapping: protein chain A and DNA chain B; position i is PDB/DeepPBS first-strand index i; no register or orientation search
- Coverage: every PDB first-strand position and all three substitutions, independent of assay coverage

## Label-free prediction summary

| protein | sequence | positions | rows | minimum predicted effect | maximum predicted effect |
|---|---|---:|---:|---:|---:|
| DBP023 | `GCAGATCTGCACATC` | 15 | 45 | -5.324860 | 0.946926 |
| DBP056 | `CGCTATCCAGAGCG` | 14 | 42 | -6.065019 | 3.282403 |
| DBP062 | `CGCGATGCTTCTCG` | 14 | 42 | -6.221465 | 3.693817 |

All three structures passed chain, duplex, output-shape, probability, and first-strand sequence checks. No experimental label was used for preprocessing, mapping, sign, checkpoint, or failure handling. Evaluation is prohibited until this artifact and report are committed.
