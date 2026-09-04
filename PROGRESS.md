# Progress

## Completed

### v0.1 historical prototype

- Initialized standalone git repository and project skeleton.
- Created and used Python 3.13 virtual environment.
- Collected 16 real PDB protein-DNA pairs with provenance.
- Validated DNA/protein QC and preserved raw JSON under `data/raw/`.
- Built `single_mutants.csv` with 843 rows and `double_mutants.csv` with 23,904 rows.
- Built `random_negatives.csv` with 32,000 rows.
- Built `benchmark_v0_1.csv` with 56,763 rows.
- Generated figures 1 to 6, including GC distribution, mutation landscapes, proxy distribution, and a chr22 genome-scan demo.
- Generated `results/WEEKLY_PROGRESS.md`.
- Published the repository to GitHub as a public remote.

### v0.2 curated benchmark

- Completed benchmark v0.2 scientific audit and documented v0.1 assumptions in `docs/SCIENTIFIC_AUDIT.md`.
- Curated all 16 historical PDB pairs in `metadata/pdb_pair_curation.csv`.
- Corrected PDB structural evidence versus quantitative specificity ground-truth annotations.
- Built `dbp_target_pairs_v0_2.csv` with 8 strict core benchmark pairs.
- Rebuilt `benchmark_v0_2.csv` with 25,761 rows: 393 single mutants, 9,360 double mutants, and 16,000 random negatives.
- Generated v0.2 figures under `results/figures/v0_2/` plus figure interpretation notes.
- Added sequence-only proxy positional-bias analysis in `analysis/05_proxy_position_bias.py`.
- Added a small Layer C experimental specificity pilot from JASPAR CORE and UniProt with 5 proteins and 1,209 PFM-derived DNA score records.
- Generated `results/BENCHMARK_QUALITY_REPORT.md` and `results/WEEKLY_PROGRESS_V0_2.md`.
- Created local tag `v0.2-benchmark-freeze` at the v0.2 benchmark state.

### v0.3 designed DBP experimental specificity benchmark

- Programmatically collected GEO GSE237017 series and GSM metadata for designed DBP uPBM experiments.
- Downloaded 12 processed 7-mer files and 12 raw spot-data files, with SHA256 and file size provenance in `metadata/v0_3/gse237017_file_manifest.csv`.
- Parsed processed uPBM files into `data/interim/gse237017/upbm_7mers_long.parquet`, explicitly expanding the reverse-complement companion 7-mer column.
- Verified complete 7-mer coverage for all 12 samples: 16,384 unique 7-mers per sample after reverse-complement expansion.
- Built `data/processed/v0_3/designed_dbp_upbm_v0_3.parquet` with 114,688 protein-7mer experimental measurements.
- Recovered all 7 designed DBP protein sequences and intended target DNA sequences from official Nature supplementary material.
- Generated replicate QC, reverse-complement QC, target-derived 7-mer rank summaries, sequence-only baseline comparisons, and v0.3 figures.
- Created v0.3 dataset card, PBM score definitions, model split plan, progress report, and GO/NO-GO report.

### v0.3.1 scientific audit and correction

- Re-read the paper/source-data definition for Extended Data Fig. 8 uPBM motif percentiles.
- Downloaded and preserved Nature Source Data Fig. 4 and Source Data Extended Data Fig. 8 under `data/raw/v0_3_1/`.
- Reproduced published uPBM motif percentiles for all 7 designed DBPs within a 2 percentile-point tolerance; maximum error was 1.2535 points.
- Separated original design target, experimental assay target, and PBM evaluation motif in `metadata/v0_3_1/designed_dbp_target_definitions.csv`; DBP48 is now original target I but assay/PBM reference C.
- Added explicit reverse-complement canonicalization and rebuilt v0.3.1 oriented and RC-class benchmark files.
- Confirmed 114,688 oriented rows collapse to 57,344 protein-RC-class experimental units.
- Recomputed RC-aware sequence-only baseline against paper motifs.
- Corrected the old "140 failure cases" interpretation: 140 is now documented as top examples only, while the v0.3.1 total sequence-vs-experiment disagreement count is 1,515.
- Added empirical replicate agreement/noise ceiling table and figure.
- Added designed DBP protein sequence clusters and target/motif groups for future split control.
- Added v0.3.1 portability/schema/reproduction tests; `pytest` reports 29 passed.
- Generated `results/v0_3_1/V0_3_1_VALIDATION_REPORT.md` with final gate `GO TO V0.4`.

### v0.4 strong protein-conditioned baseline arena

- Created local freeze tag `v0.3.1-benchmark-freeze` before starting v0.4 work.
- Added official external repositories as tracked dependencies for DeepPBS, NA-MPNN, and the designed-DBP structure source repository.
- Wrote baseline feasibility audit, task harmonization notes, failure-resolution definition, and natural PBM/uPBM control plan under `docs/v0_4/`.
- Built `metadata/v0_4/designed_dbp_structure_manifest.csv`: public checked structures are available for DBP35 and DBP48 only; DBP1, DBP3, DBP5, DBP6, and DBP9 are not structure-evaluable in this pass.
- Ran official NA-MPNN specificity inference on the bundled example, DBP35 theoretical complex model, and DBP48/8TAC experimental structure.
- Detected that `8tac` appears in NA-MPNN split files, so DBP48/8TAC is recorded as diagnostic and not zero-shot.
- Recorded DeepPBS as not fairly evaluable in the current Windows environment because the official preprocessing stack requires additional Linux-oriented structure-processing dependencies.
- Implemented `SimpleProteinConditionalBaseline` as an untrained Tier 1 protein-conditioned scaffold; no model was trained on GSE237017.
- Evaluated all baselines as per-protein RC-class ranking tasks and wrote v0.4 tables under `results/v0_4/tables/`.
- Generated six v0.4 figures under `results/v0_4/figures/`.
- Generated `results/v0_4/BASELINE_GAP_ANALYSIS.md`, `results/v0_4/NEW_MODEL_GO_NO_GO.md`, and `results/v0_4/V0_4_PROGRESS.md`.
- Current v0.4 gate: `CONDITIONAL GO`.

### v0.4.1 natural PBM training benchmark and DeepPBS Linux audit

- Created local freeze tag `v0.4-baseline-diagnostic-freeze` before starting v0.4.1 work.
- Audited natural PBM/uPBM sources and selected UniPROBE processed contiguous 8-mer E-score profiles as the primary natural training source.
- Downloaded 11 UniPROBE publication-level 8-mer archives with SHA256 provenance in `metadata/v0_4_1/natural_pbm_files.csv`.
- Parsed 112 profile files into `data/interim/v0_4_1/natural_pbm_long.parquet`.
- Recovered conservative full-length UniProt reference sequences for 57 train-ready natural protein/construct IDs; construct-level sequence matching remains unresolved and is marked `sequence_match_to_assay=false`.
- Built `data/processed/v0_4_1/natural_pbm_benchmark_v0_4_1.parquet` with 57 proteins and 1,875,072 protein-RC-class 8-mer units.
- Generated 40% protein-cluster-aware split: train 39, validation 9, natural_test 9 proteins.
- Trained a lightweight `SimpleProteinConditionalBaseline_composition_ridge` only on natural PBM train proteins; this is not the proposed method.
- SimplePC macro median Spearman: natural_test 0.301 and designed external GSE237017 0.362.
- SimplePC resolved 333/1,515 pre-registered v0.3.1 sequence-vs-experiment disagreement candidates.
- Added DeepPBS Docker/Linux wrapper and provenance, but official DeepPBS execution was not run because the current host lacks Docker/WSL runtime.
- Added v0.4.1 regression tests; `pytest` reports 61 passed.
- Current v0.4.1 gate: `WAIT FOR STRONGER BASELINE`.

### v0.4.2 strong baseline closure

- Created a v0.4.2 construct audit for the 57 natural UniPROBE proteins.
- Confirmed that no assay-aligned construct sequences were recovered from the local provenance, so the assay-aligned natural construct benchmark remains empty.
- Preserved `FULL_LENGTH_REFERENCE` as a sensitivity benchmark only.
- Recorded DeepPBS Linux runtime status and kept the official Docker/WSL limitation explicit in provenance.
- Trained a frozen ESM-2 `esm2_t12_35M_UR50D` protein-conditioned baseline on the natural PBM train split.
- FrozenPLM macro median Spearman: natural_test 0.316 and designed_external 0.153.
- Reanalyzed the 1,515 pre-registered v0.3.1 sequence-vs-experiment disagreement candidates with a frozen resolution protocol.
- Disagreement resolution counts: sequence k-mer 309/1,515, SimplePC 333/1,515, FrozenPLM 159/1,515, NA-MPNN diagnostic 50/398 evaluable candidates, DeepPBS 0 evaluable candidates.
- Defined a 263-sequence common hard set where high experimental PBM E-score candidates remain low-ranked by all complete core baselines.
- Added designed DBP difficulty diagnostics covering motif-distance regimes, protein train-set similarity, ESM embedding distance, score distribution shape, and per-protein baseline performance.
- Wrote `results/v0_4_2/FINAL_STRONG_BASELINE_GATE.md` and `results/v0_4_2/V0_4_2_VALIDATION_REPORT.md`.
- Added v0.4.2 regression tests; `pytest` reports 72 passed.
- Current v0.4.2 gate: `WAIT - BENCHMARK STILL INCOMPLETE`.

## In progress

- Preparing a supported Linux/Docker execution path for DeepPBS official preprocessing and prediction.
- Expanding natural PBM protein sequence curation from full-length UniProt references to assay construct sequences.
- Improving structure-aware baseline coverage for all seven designed DBPs.

## Next step

- Run DeepPBS in a supported Linux/Docker/WSL environment using the wrapper in `external/deeppbs/`.
- Replace v0.4.1 proxy protein clustering with MMseqs2/CD-HIT once a Linux runtime is available.
- Curate experimental PBM construct sequences where source publications provide them.
- Keep GSE237017 designed DBPs as an external test set rather than mixing them into natural DBP training.

## Known limitations

- v0.1/v0.2 scores are sequence-only proxy metrics.
- v0.3 PBM E-scores are experimental 7-mer specificity/enrichment scores, not Kd, binding free energy, binding probability, or absolute cross-protein affinity.
- Sequence-only proxy mutation landscapes can contain positional artifacts from k-mer overlap.
- The genome scan is a chr22 candidate retrieval demo, not a full-GRCh38 screen or off-target prediction.
- PDB structural pairs do not provide quantitative specificity ground truth.
- JASPAR pilot scores are PFM-derived PWM log2-odds values, not raw PBM/HT-SELEX enrichment.
- GSE237017 processed PBM tables are 7-mer level; intended full target DNA sequences are longer than 7 bp and are summarized by overlapping 7-mers only.
- Natural-to-designed OOD evaluation may be confounded by assay shift unless natural PBM/uPBM controls are added.
- No calibrated uncertainty or protein-conditioned binding model is available yet.
- v0.4 NA-MPNN results are diagnostic for two proteins only; DBP48/8TAC is not zero-shot because of detected split overlap.
- DeepPBS has not been fairly evaluated yet.
- v0.4.1 natural PBM uses 8-mer UniPROBE E-scores, while GSE237017 designed uPBM uses 7-mer E-scores; natural-to-designed comparisons remain confounded by assay/k-mer processing differences.
- v0.4.1 natural protein sequences are full-length UniProt references, not confirmed assay constructs.
- v0.4.1 SimplePC is a low-capacity baseline, not the proposed model and not a replacement for strong structure-aware baselines.
- v0.4.2 assay-aligned natural PBM construct sequences were not recovered from the current local provenance, so the primary construct-aware benchmark remains empty.
- v0.4.2 DeepPBS still lacks a runnable Linux runtime on this host; the official example remains unexecuted locally.
