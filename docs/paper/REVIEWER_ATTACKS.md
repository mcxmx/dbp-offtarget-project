# Reviewer Attack Simulation

## Reviewer 1: ML / benchmark

**Attack:** N_protein=7 is underpowered; development exposure and heterogeneous baselines compromise fairness; candidate-level N inflates significance.

**Current answer:** protein is explicitly the independent unit, all per-protein values are reported, exposure is marked, and unmatched transfer rows are context only. Evidence: `docs/v0_8_sota/SOTA_BENCHMARK_PROTOCOL.md`, `results/v0_9_external/external_dataset_audit.csv` and frozen v0.5/v0.6 manifests.

**New experiment needed:** independent assay-matched cohort with more than seven proteins.

**Severity:** high.

## Reviewer 2: Protein-DNA specificity

**Attack:** PBM E-score is not affinity; register and strand conventions may drive results; natural and designed assays are not comparable.

**Current answer:** manuscript calls scores processed ranking measurements, not affinity; RC classes and register rules are explicit; natural/design differences are reported as a hypothesis. Evidence: `docs/v0_9_external/ASSAY_SPECIFIC_VALIDATION_RULES.md`, v0.8 protocol and replicate reference.

**New experiment needed:** exact construct/target metadata and a matched external assay.

**Severity:** medium-high.

## Reviewer 3: Protein design / structure modeling

**Attack:** DeepPBS/NA-MPNN cover only 2/7, AF3 is missing, and sequence-model limitations may be unsurprising.

**Current answer:** agree with the coverage limitation; structure rows are supplementary diagnostics, not a comprehensive ranking; no universal computational claim is made. Evidence: `docs/v0_8_sota/AF3_FEASIBILITY_AUDIT.md`, `results/v0_8_sota/sota_benchmark_master.csv`.

**New experiment needed:** preregistered, input-compatible structure benchmark across the independent cohort.

**Severity:** high, not fatal for a scoped benchmark paper.

## Cross-cutting fatal risk

Without independent confirmatory data, a reviewer can reasonably call the result a single-cohort benchmark. The appropriate response is to label the manuscript developmental and state `BLOCKED ON INDEPENDENT CONFIRMATORY DATA`, not to relabel GSE237017 as external validation.
