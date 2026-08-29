# Progress

## Completed

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
- Completed benchmark v0.2 scientific audit and documented v0.1 assumptions in `docs/SCIENTIFIC_AUDIT.md`.
- Curated all 16 historical PDB pairs in `metadata/pdb_pair_curation.csv`.
- Corrected PDB structural evidence versus quantitative specificity ground-truth annotations.
- Built `dbp_target_pairs_v0_2.csv` with 8 strict core benchmark pairs.
- Rebuilt `benchmark_v0_2.csv` with 25,761 rows: 393 single mutants, 9,360 double mutants, and 16,000 random negatives.
- Generated v0.2 figures under `results/figures/v0_2/` plus figure interpretation notes.
- Added sequence-only proxy positional-bias analysis in `analysis/05_proxy_position_bias.py`.
- Added a small Layer C experimental specificity pilot from JASPAR CORE and UniProt with 5 proteins and 1,209 PFM-derived DNA score records.
- Generated `results/BENCHMARK_QUALITY_REPORT.md` and `results/WEEKLY_PROGRESS_V0_2.md`.

## In progress

- Expanding beyond structural complexes to raw PBM / HT-SELEX / CIS-BP-style quantitative specificity datasets.
- Designing protein-conditioned model interfaces for future calibration work.

## Next step

- Replace or complement the JASPAR PFM-derived pilot with raw experimental specificity measurements where accessible.
- Add CIS-BP / PBM / HT-SELEX sources and protein-family / DNA-similarity split logic.
- Attach a protein-conditioned scoring backend when a suitable model is available.

## Known limitations

- Current scores are sequence-only proxy metrics.
- Sequence-only proxy mutation landscapes can contain positional artifacts from k-mer overlap.
- The genome scan is a chr22 candidate retrieval demo, not a full-GRCh38 screen or off-target prediction.
- PDB structural pairs do not provide quantitative specificity ground truth.
- JASPAR pilot scores are PFM-derived PWM log2-odds values, not raw PBM/HT-SELEX enrichment.
- No calibrated uncertainty or protein-conditioned binding model is available yet.
