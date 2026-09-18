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

## v1.3.2 External Generalization & Identifiability (2026-09-18)

### Completed

- Froze `V1_3_2_EXTERNAL_VALIDATION_CONTRACT.md` before reading external results. Canonical Watson-Crick mutations remain primary; mismatches remain secondary; technical spot splits are not biological replicates.
- Downloaded and cached official Afek et al. SaMBA sources: Springer MOESM4/MOESM7 and GEO `GSE156375_RAW.tar`. Provenance and one-time hashes are in `V1_3_2_EXTERNAL_DATA_PROVENANCE.md` and `artifacts/cache/file_manifest.json`.
- Parsed the calibration workbook into `data/processed/v1_3_2_samba_spots.parquet`: 22,026 spot rows, 12 source sites, 7 TFs. Ten sites pass exact reconstruction of every published median; Ets1 site 2 and site 4 remain in the raw table but are excluded from endpoints because their published medians cannot be reconstructed without an undocumented block-selection rule.
- Ran the frozen seed-1301 spot split, 100-seed secondary stability analysis, 2,000-permutation nulls, landscape variance allocation, and publication-oriented Figure 3 source/render. The canonical TF-level medians are global 0.947, position 0.881, identity residual 0.823, and identity pairwise 0.814.
- Audited official SAMPDI-3D T227 sources. The paper/supplement reports aggregate values only; no auditable row-level T227 or per-row SAMPDI/FoldX predictions were recovered. The 17-versus-18 TF discrepancy is recorded; no T227 decomposition is fabricated.
- Added `manuscript/OUTLINE.md` and external audit/validation reports. No model was trained and the v1.4 gate remains closed.

### Bugs / source issues

- A concurrent manifest write initially omitted the MOESM7 cache entry; the already-computed SHA256, file size, mtime, URL, and timestamp were restored directly into the persistent manifest without rereading the workbook.
- The calibration workbook does not provide reliable array/batch identifiers. Processed rows therefore use `array_id=unknown_array` and record `array_id_source=not_supplied_by_source_workbook`.
- Ets1 site 2/site 4 sequence joins contain published-median inconsistencies. They are explicitly excluded by source-integrity status rather than mapped by outcome.

### Stable conclusions

- SaMBA canonical identity effects are highly repeatable at the technical spot-split level (median residual Spearman 0.823; median pairwise accuracy 0.814), so identity information is not universally unobservable under this assay.
- This result is not a biological replicate noise ceiling and does not identify all PBM/competition estimand or context differences.
- T227 natural-TF generalization remains not evaluable from official row-level data. The designed benchmark therefore retains its prospective position-dominant finding but cannot yet claim natural-TF universality.

### Next single priority

Recover a legitimate row-level natural-TF mutation benchmark or an independent replicate-resolved competition assay before opening v1.4 identity-aware model development.

## v1.3.3 Natural-TF Structural Model Challenge (2026-09-18)

### Completed

- Froze `V1_3_3_NATURAL_TF_CHALLENGE_CONTRACT.md` before natural-TF model/label performance inspection.
- Audited all 12 SaMBA sites across seven TFs using official RCSB experimental structures and a deterministic label-blind DNA mapping. Four ETS1 sites mapped uniquely to PDB 2STT; eight sites were excluded for mapping/coverage ambiguity. Structure files and one-time hashes are in `data/raw/v1_3_3_external/structures/` and the persistent manifest.
- Committed the structure/mapping freeze as `2e29f3c`.
- Ran the frozen official DeepPBS preprocessing in the documented Linux CPU environment. The fixed 2STT NMR input produced no NPZ; no alternate structure, register, orientation, checkpoint, or implementation was substituted.
- Committed the unscored prediction artifact as `e44e29e`. It contains 132 eligible canonical mutation keys with explicit `NOT_REPRODUCED_DEEPPBS_PREPROCESSING_NO_NPZ` status and no experimental effects.
- After the freeze, revealed canonical SaMBA labels and wrote the natural-TF decomposition tables. All primary natural-TF DeepPBS metrics are `NOT_EVALUABLE`; no score was imputed.
- Generated `figure4_natural_vs_designed_identity.png/pdf` and traceable source data. The figure distinguishes SaMBA technical repeatability, model performance, prospective holdout, and natural-TF not-evaluable status.
- Added the SAMPDI leakage status and an explicit NOT_EVALUABLE orthogonal SaMBA validation table. T227 was not reopened.

### Stable conclusions

- SaMBA technical identity repeatability remains high under its within-assay spot split, but the natural-TF model challenge currently has no valid DeepPBS performance estimate.
- The natural-TF evidence cannot distinguish designed-protein domain shift from a general structure-model identity limitation. Structural coverage is limited and all eligible sites belong to ETS1.
- v1.4 remains closed; no GNN was trained and no frozen v1.3/v1.3.1/v1.3.2 endpoint was changed.

### Bugs / reproducibility findings

- Official DeepPBS preprocessing silently skipped the fixed 2STT NMR input while returning success. This is recorded as a reproducibility failure and preserved as NA rows rather than being treated as a negative model result.

### Next single priority

Obtain a reproducible experimental-structure-compatible natural-TF inference path or an independent row-level mutation benchmark; do not open v1.4 model development from the current NOT_EVALUABLE result.

## v1.3.4 Technical Rescue + Paper Consolidation (2026-09-18)

### Completed

- Wrote `V1_3_4_TECHNICAL_RESCUE_PROTOCOL.md` and `reports/v1.3/V1_3_4_DEEPPBS_FAILURE_AUDIT.md` before rerunning any prediction.
- Confirmed that frozen 2STT is a 25-model NMR ensemble: 2,670 atoms per model, protein chain C, complementary DNA chains A/B, and no alternate locations, waters, or heteroatom records in the audited atom blocks.
- Deterministically extracted MODEL 1 to `data/processed/v1_3_4/2STT_model1.pdb` without minimization, remodeling, mutation, fitting, or register changes. Its SHA256 was added once to the persistent manifest.
- Reused the cached DBP035 positive-control output. The one permitted local MODEL 1 rescue failed at import with `ModuleNotFoundError: torch_cluster` and produced no NPZ. Historical Docker/WSL unavailability was reused from the existing audit; no alternate route, conformer, structure, or webserver result was substituted.
- Closed the bounded natural-TF rescue with primary status `NOT_EVALUABLE_PRIMARY_UNCHANGED`. No post-reveal prediction table was fabricated.
- Generated `results/v1_3/variance_decomposition_per_protein.tsv` and `results/v1_3/position_only_counterfactual.tsv` from frozen experimental and prediction tables. These distinguish landscape variance allocation from predictive performance.
- Generated publication-oriented Figures 1-5 and the full `manuscript/MANUSCRIPT_DRAFT.md`, `manuscript/SUPPLEMENT_PLAN.md`, and `manuscript/SUBMISSION_GAP_ANALYSIS.md`.

### Stable v1.3.4 conclusions

- Experimental landscapes contain a nonzero within-position identity component, while DeepPBS often allocates a larger share of predicted variance to between-position effects.
- Position-only counterfactuals can preserve aggregate global ranking, so aggregate performance can mask loss of nucleotide-identity discrimination.
- The locked holdout continues to support `position >> identity` (global 0.085, position 0.560, identity residual 0.001, pairwise 0.548 median).
- SaMBA identity repeatability is high in a distinct technical assay (median residual Spearman 0.823; pairwise 0.814) and is not a competition noise ceiling.
- Protein conditionality and DNA mutation resolution remain orthogonal failure axes.

### Gate and stop condition

`V1.4 GATE: CLOSED.` The evidence supports manuscript consolidation, but the training set beyond seven designed proteins and replicate-resolved competition identifiability are not yet defensible for identity-aware model development. No new neural network or GNN was trained. v1.3.4 is complete; the next action is manuscript refinement/submission preparation or a separately approved v1.4 project, not another v1.3 benchmark expansion.
