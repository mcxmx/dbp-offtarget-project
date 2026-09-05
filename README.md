# DBP Off-target Prototype

Reproducible research prototype for genome-wide off-target risk assessment and uncertainty calibration for sequence-specific DNA-binding proteins.

Current status: v0.5 target-conditioned specificity model family implemented;
one fixed-fold engineering smoke training completed. The complete four-fold
primary evaluation and final proposed model have not been run.

Current scope:
- public structural protein-DNA pairs from RCSB PDB, curated as a structural cognate layer
- single- and double-mutant DNA perturbation benchmarks
- GC-matched and random negative controls
- sequence-only proxy baselines for pipeline sanity checks
- chr22 genome candidate retrieval prototype
- designed DBP uPBM experimental specificity benchmark from GEO GSE237017
- baseline feasibility and diagnostic evaluation for DeepPBS / NA-MPNN / simple protein-conditioned scaffold
- natural UniPROBE PBM 8-mer training benchmark for first protein-conditioned baseline sanity checks
- v0.5 matched model family (`M0` to `M3`) with RC-invariant scoring, deterministic pairwise ranking training, and one legal-fold smoke test

NO FINAL PROPOSED MODEL TRAINED YET.

DeepPBS status: the official Linux workflow was reproduced in the configured
Ubuntu VM. DBP35 and DBP48 are evaluable; coverage is 2/7 because no
legitimate public/project protein-DNA structure was found for the other five
designed DBPs.

Next stage: use DeepPBS as a completed partial structure-aware diagnostic and
resolve structure coverage/overlap limitations before starting a final
proposed model.

## v0.5 Target Contract

- Preserves independently sourced intended design target, experimental assay
  reference, and PBM-derived motif as separate fields.
- Uses the existing four protein-sequence clusters for the primary
  leave-one-protein-cluster-out split.
- Adds a stricter three-component split using the combined protein/target/motif
  leakage graph.
- Groups all future DNA-level partitions by canonical reverse-complement
  equivalence class.
- Prohibits random splitting of the 57,344 protein-RC-class rows.
- The v0.5 matched implementation is an engineering/smoke-test family, not the final proposed model; no complete four-fold primary evaluation has been run.

## v0.4.2 Current Results

- Natural PBM construct audit: no assay-aligned construct sequences were recovered from the current local provenance for the 57 UniPROBE proteins, so the primary assay-aligned construct benchmark remains empty and `FULL_LENGTH_REFERENCE` is kept only as a sensitivity benchmark.
- FrozenPLMProteinConditionalBaseline: frozen ESM-2 `esm2_t12_35M_UR50D` mean-pooled protein embeddings with a small ridge head; natural_test macro median Spearman 0.316 and designed_external macro median Spearman 0.153.
- Compared with the prior SimplePC baseline, FrozenPLM improves natural held-out ranking slightly but is weaker on the designed DBP external set.
- Best prior designed sequence-only baseline remains 0.232 macro median Spearman; empirical designed uPBM replicate reference remains about 0.591.
- Reanalyzed the pre-registered 1,515 v0.3.1 sequence-vs-experiment disagreement candidates; SimplePC resolves 333, sequence k-mer resolves 309, FrozenPLM resolves 159, and 263 candidates remain common high-experiment/low-core-baseline cases.
- DeepPBS official `5x6g` example passed in the Ubuntu VM. DBP35 and DBP48 also completed the official preprocessing/inference workflow and were evaluated over 8,192 RC-class 7-mers each.
- DeepPBS diagnostic macro median Spearman is 0.159 across 2/7 evaluable designed DBPs: DBP35 0.040 and DBP48 0.278. This is not a seven-protein generalization estimate.
- DeepPBS-derived scores are PWM/log-probability ranking proxies, not affinity, Kd, binding probability, or calibrated specificity scores.
- The completed artifact is `results/v0_4_2/DEEPPBS_COMPLETION_REPORT.md`; the earlier host-limited WAIT report is retained as historical output.
- v0.4.2 DeepPBS gate: `WAIT FOR STRONGER STRUCTURE-COVERED BASELINE`.

## v0.4.1 Current Results

- Natural PBM source: UniPROBE processed contiguous 8-mer E-score profiles.
- Final natural benchmark: 57 proteins, 7 coarse protein-family classes, 7 species, and 1,875,072 protein-RC-class 8-mer units.
- 40% protein-cluster split: train 39, validation 9, natural_test 9 proteins.
- SimpleProteinConditionalBaseline: lightweight composition/ridge baseline trained only on natural PBM train proteins; not the proposed method.
- SimplePC macro median Spearman: natural_test 0.301; designed external GSE237017 0.362.
- Best prior designed sequence-only baseline: 0.232 macro median Spearman.
- Designed uPBM empirical replicate Spearman reference: about 0.591.
- DeepPBS: official Linux runtime and the two legally evaluable designed-DBP inputs were completed in the Ubuntu VM; coverage remains 2/7.
- v0.4.1 gate: `WAIT FOR STRONGER BASELINE`.

## v0.4.1 Run Order

1. `.\.venv313\Scripts\python analysis/v0_4_1/00_natural_pbm_source_audit_and_download.py`
2. `.\.venv313\Scripts\python analysis/v0_4_1/01_build_natural_pbm_benchmark.py`
3. `.\.venv313\Scripts\python analysis/v0_4_1/02_recover_natural_pbm_sequences.py`
4. `.\.venv313\Scripts\python analysis/v0_4_1/03_finalize_natural_pbm_benchmark.py`
5. `.\.venv313\Scripts\python analysis/v0_4_1/04_train_simple_pc_baseline.py`
6. `.\.venv313\Scripts\python analysis/v0_4_1/05_deeppbs_linux_reproduction_audit.py`
7. `.\.venv313\Scripts\python analysis/v0_4_1/06_v0_4_1_reports_and_figures.py`
8. `.\.venv313\Scripts\python -m pytest`

## v0.4 Current Results

- Fixed benchmark: v0.3.1 GSE237017 designed-DBP uPBM, 7 proteins and 57,344 protein-RC-class experimental units.
- Tier 0 sequence-only best baseline: `sequence_kmer3`, macro median Spearman 0.232.
- Empirical replicate agreement reference: median E-score replicate Spearman 0.591.
- DeepPBS: official repository checked, not fairly evaluable in the current Windows environment.
- NA-MPNN: official specificity inference ran for DBP35 and DBP48 only; DBP48/8TAC has split-overlap risk and is diagnostic, not zero-shot.
- SimpleProteinConditionalBaseline: implemented as an untrained Tier 1 interface scaffold; not reported as a model result.
- v0.4 gate: `CONDITIONAL GO`, because sequence-only gap is clear but strong structure-aware baseline coverage is incomplete.

## v0.4 Run Order

1. `.\.venv313\Scripts\python analysis/v0_4/00_baseline_feasibility_and_manifests.py`
2. Run official NA-MPNN specificity inference for evaluable structures, preserving outputs under `results/v0_4/external_runs/`.
3. `.\.venv313\Scripts\python analysis/v0_4/01_build_v0_4_predictions_and_evaluation.py`
4. `.\.venv313\Scripts\python analysis/v0_4/02_generate_v0_4_figures.py`
5. `.\.venv313\Scripts\python analysis/v0_4/03_write_v0_4_reports.py`
6. `.\.venv313\Scripts\python -m pytest -q`

## v0.3.1 Validation Run Order

1. `.\.venv313\Scripts\python analysis/v0_3_1/00_prepare_target_definitions_and_manifest.py`
2. `.\.venv313\Scripts\python analysis/v0_3_1/01_reproduce_paper_upbm_percentiles.py`
3. `.\.venv313\Scripts\python src/build_designed_dbp_benchmark_v0_3_1.py`
4. `.\.venv313\Scripts\python analysis/v0_3_1/02_rc_aware_sequence_baseline.py`
5. `.\.venv313\Scripts\python analysis/v0_3_1/03_disagreement_and_noise_ceiling.py`
6. `.\.venv313\Scripts\python analysis/v0_3_1/04_protein_target_grouping.py`
7. `.\.venv313\Scripts\python -m pytest -q`
8. `.\.venv313\Scripts\python analysis/v0_3_1/05_write_validation_report.py`

## v0.3.1 Key Numbers

- Designed DBPs: 7
- uPBM samples: 12
- Oriented 7-mers per protein: 16384
- RC equivalence classes per protein: 8192
- Total oriented protein-7mer rows: 114688
- Total protein-RC-class experimental units: 57344
- Extended Data Fig. 8 percentile reproduction: PASS for 7/7 proteins
- Tests: 29 passed

## v0.3 Run Order

1. `.\.venv313\Scripts\python src/collect_gse237017.py`
2. `.\.venv313\Scripts\python src/parse_upbm.py`
3. `.\.venv313\Scripts\python analysis/v0_3/01_upbm_qc.py`
4. `.\.venv313\Scripts\python src/collect_designed_dbp_metadata.py`
5. `.\.venv313\Scripts\python src/build_designed_dbp_benchmark_v0_3.py`
6. `.\.venv313\Scripts\python analysis/v0_3/02_experimental_figures_and_baseline.py`
7. `.\.venv313\Scripts\python analysis/v0_3/03_write_v0_3_reports.py`

## v0.2 Run order

1. `.\.venv313\Scripts\python src/curate_pdb_pairs.py`
2. `.\.venv313\Scripts\python src/build_benchmark_v0_2.py`
3. `.\.venv313\Scripts\python analysis/05_proxy_position_bias.py`
4. `.\.venv313\Scripts\python analysis/06_v0_2_figures.py`
5. `.\.venv313\Scripts\python src/collect_experimental_specificity.py`
6. `.\.venv313\Scripts\python analysis/08_experimental_specificity_proxy_baseline.py`
7. `.\.venv313\Scripts\python analysis/07_benchmark_quality_report.py`

## v0.1 Historical Run Order

1. `.\.venv313\Scripts\python src/collect_data.py`
2. `.\.venv313\Scripts\python src/validate_sequences.py`
3. `.\.venv313\Scripts\python src/preprocess.py`
4. `.\.venv313\Scripts\python src/generate_mutants.py`
5. `.\.venv313\Scripts\python src/generate_gc_matched.py`
6. `.\.venv313\Scripts\python src/preprocess.py`
7. `.\.venv313\Scripts\python src/sequence_baselines.py`
8. `.\.venv313\Scripts\python analysis/01_dataset_summary.py`
9. `.\.venv313\Scripts\python analysis/02_mutation_landscape.py`
10. `.\.venv313\Scripts\python analysis/03_sequence_baseline.py`
11. `.\.venv313\Scripts\python analysis/04_genome_scan_demo.py`

## Outputs

- `data/raw/` raw API and genome downloads
- `data/interim/` staged tables
- `data/processed/` curated tables
- `results/figures/` figures
- `results/tables/` summaries

## Notes

PDB structures are treated as structural cognate evidence, not quantitative specificity ground truth. v0.1/v0.2 sequence scores are sequence-only proxy metrics for sanity checking. They are not protein-conditioned binding predictions, binding affinities, calibrated risks, or biological specificity landscapes.

v0.3 adds real in vitro uPBM experimental 7-mer specificity measurements for designed DBPs. The primary score is PBM E-score for per-protein ranking. It must not be interpreted as Kd, binding free energy, binding probability, in vivo genomic binding, or absolute cross-protein affinity.

v0.3.1 supersedes the v0.3 target-rank interpretation by separating original design target, experimental assay target, and PBM motif reference. It also makes reverse-complement equivalence explicit and reports 57344 protein-RC-class experimental units.
