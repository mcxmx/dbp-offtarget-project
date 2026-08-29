# Scientific Audit: Benchmark v0.1 to v0.2

Retrieval/audit date: 2026-08-29

This document records the scientific audit of the existing v0.1 prototype before
changing the benchmark. The purpose of v0.2 is to make the dataset interpretable
as a small research benchmark, not to expand features.

## Current v0.1 State

- The repository contains 16 RCSB PDB-derived protein-DNA pairs in
  `data/processed/dbp_target_pairs.csv`.
- v0.1 generated target, single-mutant, double-mutant, GC-matched random,
  random DNA, sequence-only proxy scores, preliminary figures, and a chr22
  genome candidate retrieval demo.
- The code preserved basic provenance such as PDB ID, chain IDs, source URL,
  and retrieval date.

## Critical Scientific Issues Found

### PDB structure was conflated with specificity ground truth

`src/collect_data.py` set `has_specificity_ground_truth=True` for PDB-derived
records. This is scientifically incorrect. A protein-DNA complex structure
supports the existence of a structural cognate complex, but it does not by
itself provide quantitative specificity ground truth.

v0.2 must split this into separate evidence fields:

- `has_structural_cognate`
- `has_direct_dna_binding_evidence`
- `has_sequence_specificity_evidence`
- `has_quantitative_specificity_ground_truth`

Unless quantitative experimental specificity data are attached, PDB-only pairs
must have `has_quantitative_specificity_ground_truth=False`.

### Longest-chain selection is not a valid biological rule

`DECISIONS.md` and `src/collect_data.py` document/use a longest protein entity
and longest DNA entity selection heuristic. This can select the wrong chain in
multi-protein complexes or systems with multiple DNA roles.

v0.2 must curate chains using available annotation, biological mechanism, and
where feasible chain-level protein-DNA contacts. If the direct recognition chain
cannot be determined, the record should be marked `uncertain` rather than
guessed.

### Longest DNA chain is not necessarily the target DNA

Some complexes can include guide DNA, substrate DNA, damaged DNA, primers,
non-specific DNA duplexes, or engineered crystallization constructs. v0.2 must
annotate `dna_role`, with values such as `target`, `cognate_site`, `guide`,
`substrate`, `damaged_substrate`, `non_specific`, and `unknown`.

### Mechanistically distinct systems are mixed in v0.1

The current 16 PDB pairs include likely examples of sequence-specific
transcription factors, non-specific chromosomal DNA-binding proteins,
guide-dependent systems, lesion-recognition systems, transposases/nucleases,
multi-protein complexes, and designed binders. These should not be pooled as one
homogeneous specificity benchmark.

v0.2 must classify records with `sequence_specificity_class` and
`recommended_use` so that unsuitable samples remain archived but are excluded
from the core specificity benchmark.

### Sequence-only figures may be misread as biological specificity

v0.1 figures and summaries use sequence-only proxy metrics. These are useful for
pipeline sanity checks, but they are not protein-conditioned binding predictions
and not biological specificity landscapes.

v0.2 figure titles, axis labels, README text, and progress reports must use
phrases such as:

- sequence-only proxy
- sequence similarity baseline
- pipeline sanity check
- genome candidate retrieval prototype

They must not imply binding affinity, binding probability, calibrated risk, or
true off-target prediction.

### k-mer metrics can create positional artifacts

The current combined proxy includes k-mer Jaccard terms. For a fixed-length DNA
target, mutations near the middle can disrupt more overlapping k-mers than
mutations near the ends. This can create a U-shaped positional pattern even when
there is no biological specificity signal.

v0.2 must explicitly analyze this positional artifact and label the result as a
metric property, not a protein-DNA recognition result.

## v0.2 Correction Targets

- Keep all v0.1 files intact.
- Create `metadata/pdb_pair_curation.csv` with per-pair biological curation.
- Create `data/processed/dbp_target_pairs_all_curated.csv` for all historical
  PDB pairs with corrected evidence annotations.
- Create `data/processed/dbp_target_pairs_v0_2.csv` for the curated subset used
  in the v0.2 synthetic perturbation benchmark.
- Rebuild v0.2 mutant/control/scored benchmark files with versioned names.
- Replot v0.2 figures under `results/figures/v0_2/` with figure notes.
- Add an experimental specificity data pipeline as a separate data layer.

## v0.2 Audit Results

### Curation outcome

The 16 historical PDB pairs were retained in
`data/processed/dbp_target_pairs_all_curated.csv`. Only 8 records were retained
for `data/processed/dbp_target_pairs_v0_2.csv` under the strict
`core_benchmark` rule:

- 7 natural sequence-specific structural cases
- 1 designed sequence-specific structural case

Guide-dependent, lesion-specific, non-specific, and transposase/substrate cases
were preserved in the curated all-pairs table but excluded from the v0.2 core
synthetic perturbation benchmark.

### Ground-truth correction

All PDB-derived v0.2 records set:

- `has_structural_cognate=True`
- `has_quantitative_specificity_ground_truth=False`

This fixes the v0.1 mistake of treating PDB complex existence as specificity
ground truth.

### Positional proxy artifact

`analysis/05_proxy_position_bias.py` compared single-mutation position effects
for GAL4/1D66 using hamming similarity, edit similarity, 3-mer Jaccard, 4-mer
Jaccard, and the combined sequence-only proxy. Hamming and edit similarity were
flat across positions, as expected for one substitution in a fixed-length
sequence. k-mer Jaccard metrics varied by position because internal bases are
covered by more overlapping k-mers than edge bases.

Observed metric ranges for this target:

- hamming similarity: 0.000000
- edit similarity: 0.000000
- 3-mer Jaccard: 0.292857
- 4-mer Jaccard: 0.358553
- combined sequence-only proxy: 0.134764

Therefore, position-dependent patterns in sequence-proxy mutation landscapes can
be metric artifacts. They must not be described as biological specificity
landscapes without protein-conditioned scoring or experimental ground truth.

