# DBP Off-target Prototype

Reproducible research prototype for genome-wide off-target risk assessment and uncertainty calibration for sequence-specific DNA-binding proteins.

Current scope:
- public structural protein-DNA pairs from RCSB PDB, curated as a structural cognate layer
- single- and double-mutant DNA perturbation benchmarks
- GC-matched and random negative controls
- sequence-only proxy baselines for pipeline sanity checks
- chr22 genome candidate retrieval prototype
- designed DBP uPBM experimental specificity benchmark from GEO GSE237017

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
