# DBP Off-target Prototype

Reproducible research prototype for genome-wide off-target risk assessment and uncertainty calibration for sequence-specific DNA-binding proteins.

Current scope:
- public structural protein-DNA pairs from RCSB PDB
- single- and double-mutant DNA benchmarks
- GC-matched and random negative controls
- sequence-only proxy baseline
- chr22 genome-scan demo

## Run order

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

All scores in this prototype are sequence-only proxy metrics. They are not protein-conditioned binding predictions.

