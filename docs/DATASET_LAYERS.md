# Dataset Layers

Retrieval/audit date: 2026-08-29

The benchmark is separated into four layers so that structural evidence,
synthetic perturbations, experimental specificity measurements, and genome
retrieval are not interpreted as the same kind of evidence.

## Layer A: Structural Cognate Dataset

Files:

- `data/processed/dbp_target_pairs_all_curated.csv`
- `metadata/pdb_pair_curation.csv`

Purpose:

- Preserve real PDB-derived protein-DNA structural pairs with provenance.
- Record PDB ID, chain IDs, UniProt ID when available, paper DOI/PMID, source
  URL, retrieval date, biological mechanism, DNA role, and curation confidence.

What it supports:

- A structural protein-DNA complex exists for the curated pair.
- Chain-level contact evidence can support direct DNA binding.

What it does not support by itself:

- Quantitative specificity ground truth.
- Binding affinity, calibrated risk, or genome-wide off-target prediction.

## Layer B: Synthetic Perturbation Dataset

Files:

- `data/processed/single_mutants_v0_2.csv`
- `data/processed/double_mutants_v0_2.csv`
- `data/processed/random_negatives_v0_2.csv`
- `data/processed/benchmark_v0_2.csv`
- `data/processed/benchmark_v0_2_scored.csv`

Purpose:

- Stress-test the pipeline around curated target DNA sequences.
- Enumerate single-base and double-base substitutions.
- Generate length-matched GC-controlled and fully random negative controls.

What it supports:

- Reproducible sequence perturbation analysis.
- Sequence-only proxy sanity checks.

What it does not support by itself:

- Biological specificity landscapes.
- Protein-conditioned binding predictions.
- Experimental calibration metrics.

## Layer C: Experimental Specificity Dataset

Files:

- `data/processed/experimental_specificity_small.csv`
- `data/processed/v0_3/designed_dbp_upbm_v0_3.parquet`
- `metadata/v0_3/gse237017_samples.csv`
- `metadata/v0_3/gse237017_file_manifest.csv`
- `metadata/v0_3/designed_dbp_sequences.csv`
- `metadata/v0_3/designed_dbp_targets.csv`

Purpose:

- Store experimentally measured DNA preference records from PBM, HT-SELEX,
  CIS-BP/JASPAR-like curated resources, Kd/competition assays, or related
  quantitative sources.
- Store GSE237017 designed DBP uPBM 7-mer specificity measurements as the first
  raw experimental designed-DBP external benchmark.

Required columns:

- `protein_id`
- `protein_name`
- `protein_sequence`
- `dna_sequence`
- `experimental_score`
- `score_type`
- `experiment_type`
- `source_database`
- `source_id`
- `paper_doi`
- `source_url`
- `notes`

What it supports:

- Quantitative validation once score meaning and assay type are clear.
- Per-protein ranking evaluation for designed DBP 7-mer PBM specificity.
- External/OOD evaluation of protein-conditioned models trained on natural DBPs.

What it does not support automatically:

- Cross-assay score normalization.
- Family-level generalization claims without leakage-controlled splits.
- Cross-protein absolute affinity comparison.
- Direct full-target affinity when only overlapping 7-mer scores are available.

v0.3 source:

- GEO GSE237017 is Layer C experimental specificity ground truth for designed
  DBPs. Its primary score is `PBM E-score`, stored as
  `experimental_score_primary` after within-protein replicate consensus.

Important distinction:

- The v0.2 JASPAR pilot is PFM-derived motif scoring and is not equivalent to
  GSE237017 raw experimental PBM measurements. It remains useful for code-path
  testing, but it should not be mixed with v0.3 uPBM records as if the score
  scales were the same.

## Layer D: Genome Candidate Dataset

Files:

- `results/tables/genome_candidates_demo.csv`

Purpose:

- Retrieve near-match genomic candidate loci for selected target DNA sequences.

Current scope:

- GRCh38 chr22 demo only.
- Seed-based retrieval with top-candidate truncation.

Correct interpretation:

- Genome candidate retrieval prototype.

Incorrect interpretation:

- Not off-target prediction.
- Not exhaustive GRCh38 screening.
- Not protein-conditioned risk assessment.
