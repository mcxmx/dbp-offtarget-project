# v1.3.3 Natural-TF Prediction Freeze

Status: **STAGE A FROZEN**, 2026-09-18.

## Frozen artifact

- Prediction file: `results/v1_3/natural_tf_predictions_unscored.tsv`
- Rows: 132 canonical mutation keys (four ETS1 sites, 11 positions, three
  non-wild-type bases per position).
- Model: official DeepPBS checkout at commit
  `8bfb211dd67f02877841f6f33aa493ddf7daedf9`, bundled ensemble checkpoints
  `828, 529, 173, 898, 820`.
- Structure/mapping freeze commit: `2e29f3c`.
- Input mapping: `results/v1_3/natural_tf_structure_mapping.tsv`.
- Prediction-stage input: frozen 2STT experimental NMR structure and the
  label-blind SaMBA mutation keys. No SaMBA effect, rank, heatmap, or correlation
  was loaded.

## Inference status

The official DeepPBS preprocessing returned without producing an NPZ for the
fixed 2STT NMR input in the documented Linux CPU environment. The unscored
file therefore preserves every eligible mutation key with
`predicted_effect=NA` and
`prediction_status=NOT_REPRODUCED_DEEPPBS_PREPROCESSING_NO_NPZ`. This is a
reproducibility failure, not a performance result. No alternate PDB, structure
register, orientation, checkpoint, or DeepPBS implementation was substituted.

## Reveal boundary

This file was written before joining SaMBA effects. Stage B may report only
`NOT_EVALUABLE` for the primary natural-TF model challenge unless a future
secondary method is explicitly labelled post-hoc. The primary eligibility and
metric definitions remain unchanged.

The artifact is recorded once in `artifacts/cache/file_manifest.json`; later
checks must use the existing path/size/mtime cache.
