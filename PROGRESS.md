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

## In progress

- Expanding beyond structural complexes to datasets with quantitative specificity ground truth.
- Designing protein-conditioned model interfaces for future calibration work.

## Next step

- Add CIS-BP / PBM / HT-SELEX sources and protein-family / DNA-similarity split logic.
- Attach a protein-conditioned scoring backend when a suitable model is available.

## Known limitations

- Current scores are sequence-only proxy metrics.
- The genome scan is a chr22 demo, not a full-GRCh38 screen.
- No calibrated uncertainty or true binding specificity ground truth is available yet.
