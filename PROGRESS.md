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

## In progress

- Preparing v0.4 protein-conditioned baseline experiments using v0.3 as an external designed-DBP benchmark.
- Expanding beyond GSE237017 to additional raw PBM / HT-SELEX / CIS-BP-style quantitative specificity datasets.
- Designing protein-conditioned model interfaces for future calibration work.

## Next step

- Keep GSE237017 designed DBPs as an external/OOD test set rather than mixing them into natural DBP training.
- Add CIS-BP / PBM / HT-SELEX sources and protein-family / DNA-similarity split logic.
- Attach a protein-conditioned scoring backend when a suitable model is available.

## Known limitations

- v0.1/v0.2 scores are sequence-only proxy metrics.
- v0.3 PBM E-scores are experimental 7-mer specificity/enrichment scores, not Kd, binding free energy, binding probability, or absolute cross-protein affinity.
- Sequence-only proxy mutation landscapes can contain positional artifacts from k-mer overlap.
- The genome scan is a chr22 candidate retrieval demo, not a full-GRCh38 screen or off-target prediction.
- PDB structural pairs do not provide quantitative specificity ground truth.
- JASPAR pilot scores are PFM-derived PWM log2-odds values, not raw PBM/HT-SELEX enrichment.
- GSE237017 processed PBM tables are 7-mer level; intended full target DNA sequences are longer than 7 bp and are summarized by overlapping 7-mers only.
- No calibrated uncertainty or protein-conditioned binding model is available yet.
