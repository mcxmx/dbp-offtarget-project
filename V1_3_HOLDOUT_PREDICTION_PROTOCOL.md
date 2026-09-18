# V1.3.1 Holdout Prediction Protocol

Status: **FROZEN BEFORE INFERENCE, 2026-09-17**. This protocol governs DBP023, DBP056, and DBP062, which remain `LOCKED_HOLDOUT`. Competition labels, heatmaps, ranks, position summaries, and performance may not be read during Stage 1.

## Primary method

1. Use the existing official DeepPBS checkout at commit `8bfb211dd67f02877841f6f33aa493ddf7daedf9` and the repository's established Linux CPU environment.
2. Preprocess with the existing `run/process/process_config.json` and `process_co_crystal.py --no_pwm`, exactly as used by the repository's `process_and_predict.sh` entry point.
3. Predict with the existing `run/process/pred_configs/pred_config_deeppbs.json` and unmodified upstream `run/predict.py`. Its checked-in `run/plot_scripts/txts/DeepPBS.txt` fixes the five bundled checkpoints `828`, `529`, `173`, `898`, and `820`; outputs are ensemble-averaged and strand-symmetrized by upstream code.
4. No checkpoint, model version, config, or sign may be selected using holdout labels. A different DeepPBS version would be secondary/post-hoc and cannot replace this primary prediction.

## Structure and mapping

1. Inputs are the existing design models `data/raw/v1_0_structure/design_pdbs/design_pdbs/DBP023.pdb`, `DBP056.pdb`, and `DBP062.pdb` without coordinate edits.
2. The sole amino-acid chain is protein chain A. DNA chain B contains a forward strand followed by its exact reverse complement.
3. The forward strand is the first half of chain B in PDB residue order. Its PDB order and orientation are retained exactly; reverse complement, register search, and performance-based orientation selection are forbidden.
4. Position `i` maps deterministically to DeepPBS output row `i` and PDB first-strand residue `i`. The emitted DeepPBS one-hot sequence must exactly equal the PDB first strand or the protein is marked failed.
5. All structure-defined first-strand positions are scored, irrespective of later assay coverage. Stage 2 may inner-join only exact `(protein, position, wt_base, mutant_base)` keys; unmatched predictions remain documented and cannot affect mapping.

## Mutation score and preprocessing

For every position and each of the three non-wild-type bases, the frozen score is

`predicted_effect = log(P[position, mutant] + 1e-12) - log(P[position, wild_type] + 1e-12)`.

This is the v1.1 native-PWM delta-logP definition with its original sign. No centering, calibration, normalization against labels, or protein-specific transformation is applied before freeze. Stage 2 uses the unchanged v1.3 evaluation contract for global, position `mean_abs_effect`, within-position residual, pairwise, and permutation metrics.

## Failure handling

Preprocessing or prediction failure, invalid probabilities, unexpected tensor shape, ambiguous chains, non-complementary duplex, or sequence mismatch produces `NOT_EVALUABLE` for that protein. No alternate chain, trimmed structure, register, orientation, config, checkpoint, or implementation may be tried as a function of holdout performance. Environment-only repair is allowed only if it preserves the method and this mapping.

## Two-stage boundary

Stage 1 writes only `results/v1_3/holdout_predictions_unscored.tsv` and `reports/v1.3/HOLDOUT_PREDICTION_FREEZE.md`. It must not open competition labels or calculate accuracy. The prediction file is hashed once through the persistent size/mtime-aware manifest and committed. Only after that commit may Stage 2 join the frozen competition table and create holdout metrics or experimental-versus-prediction figures.
