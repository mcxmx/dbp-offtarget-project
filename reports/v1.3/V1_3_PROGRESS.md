# v1.3 progress

## Completed

- Frozen the v1.3 data audit and evaluation contract before benchmark results. The primary endpoints were not changed in v1.3.1.
- Downloaded the official Extended Data Fig. 3 and Fig. 9 workbooks from Springer Nature, retained them under `data/raw/v1_3_external/`, and cached one SHA256 per source in the persistent file manifest.
- Parsed all ten official competition landscapes. The 414 published means exactly reconstruct the repository workbook (maximum and median absolute error both 0).
- Audited the workbooks for replicate-level values. Both provide one published summary column for the relevant panels and no replicate-valued columns; no replicate values were inferred or digitized.
- Recovered the 42-mutation, 14-position DBP35opt landscape and assay metadata from Extended Data Fig. 9e. The experimental decomposition was activated as a data-availability amendment.
- Froze the DeepPBS protocol, ran label-blind prospective inference for DBP023/056/062, wrote `results/v1_3/holdout_predictions_unscored.tsv`, and committed the Stage 1 artifact as `ff68394` before joining labels.
- After the freeze commit, evaluated the three proteins with the unchanged v1.3 contract and generated the revised Figure 2, reproducibility Figure 3, and DBP35/DBP35opt Figure 4 with traceable TSV source data.

## Inputs and outputs

- Official sources: `data/raw/v1_3_external/41594_2025_1669_MOESM16_ESM.xls` and `41594_2025_1669_MOESM21_ESM.xls`.
- Replicate audit: `results/v1_3/competition_mean_reconstruction.tsv`, `competition_replicate_availability.tsv`, `replicate_noise_ceiling.tsv`, and `assay_reproducibility_matrix.tsv`.
- Prospective prediction freeze: `V1_3_HOLDOUT_PREDICTION_PROTOCOL.md`, `results/v1_3/holdout_predictions_unscored.tsv`, and `reports/v1.3/HOLDOUT_PREDICTION_FREEZE.md`.
- Reveal outputs: `results/v1_3/holdout_decomposition.tsv`, `holdout_decomposition_summary.tsv`, `holdout_permutation_nulls.tsv`, and `holdout_predictions_revealed.tsv`.
- Figure source data: `results/v1_3/figure2_revised_data.tsv`, `figure3_reproducibility_data.tsv`, and `figure4_dbp35opt_data.tsv`.

## Main numbers

- Locked holdout DeepPBS median Spearman: global 0.085, position 0.560, and within-position residual identity 0.001 across three proteins.
- Locked holdout median within-position pairwise accuracy: 0.548, versus the fixed chance baseline of 0.5.
- Per-protein position/identity residual Spearman: DBP023 0.560/0.006; DBP056 0.455/0.001; DBP062 0.705/-0.038.
- DBP35 versus DBP35opt experimental conservation: global Spearman 0.736, position-sensitivity Spearman 0.776, and identity-residual Spearman 0.184.
- DBP35opt mean absolute landscape shift is 0.357 normalized PE/FITC units. This is context-confounded because the two assays used different target and competitor concentrations.

## Bugs found

- The Stage 2 report writer initially accessed a pandas `median` method rather than the `median` column. This affected only report rendering; the stored metrics and permutation outputs were complete. The column access was corrected and the completed evaluation was reused without rerunning permutations.
- Figure 4 initially derived display limits with NaN-sensitive extrema because wild-type cells are intentionally missing. The plotting code now uses NaN-aware extrema; source values and metrics were unchanged.

## Stable conclusions

- The three prospective holdouts directionally reproduce the development finding: DeepPBS position sensitivity is substantially stronger than within-position nucleotide identity, while pairwise identity accuracy remains close to chance.
- The holdout sample is only three proteins; this is prospective directional validation, not a broad significance claim.
- The official downloadable competition source data do not support a replicate noise-ceiling calculation. Weak identity prediction therefore cannot yet be assigned uniquely to model failure rather than assay-level identifiability.
- DBP35opt is experimentally recoverable, but the assay-context difference prevents interpreting DBP35-to-DBP35opt differences as an isolated causal protein-mutation effect. No DBP35opt model prediction exists.

## Not yet concluded

- Competition mutation-, position-, and identity-level replicate reproducibility remain `NOT_EVALUABLE_SOURCE_DATA_MEAN_ONLY`.
- The predicted DBP35-to-DBP35opt landscape-shift endpoint remains `NOT_EVALUABLE_NO_DBP35OPT_MODEL_PREDICTION`.
- The measurement gate for an identity-aware v1.4 model is not met; no new GNN was trained.

## Next single priority

Obtain the underlying competition replicate-level records from an author or public repository, or establish an independent replicate-resolved mutation benchmark, before deciding whether identity-aware model development is scientifically justified.
