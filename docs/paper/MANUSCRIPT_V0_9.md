# Sequence-Based Specificity Prediction for De Novo DNA-Binding Proteins: A Generalization-Gap Benchmark

## Abstract

Predicting DNA-binding specificity for de novo proteins requires transfer to protein sequences that were not represented during training. We establish a leakage-controlled benchmark using seven designed DNA-binding proteins and separate ranking performance from evidence that predictions change with protein or target conditioning. Existing sequence-based approaches show a reproducible benchmark-level gap under strict unseen-protein evaluation. Repairing a reverse-complement representation collision improves DNA ranking, but protein- and target-shuffle diagnostics remain effectively collapsed. Increasing DNA observations, dense supervision and local residue attention do not rescue this failure. Partial DeepPBS and NA-MPNN diagnostics are structure-conditioned but cover only two of seven exposed proteins, so they cannot support a universal claim about computational methods. Hard-case analysis identifies 871 of 1,515 cases unresolved by all current sequence models. An independent validation cohort is preregistered but pending.

## Introduction

Sequence-specific DNA recognition is a central requirement for programmable protein design. De novo DNA-binding proteins differ from natural transcription factors in sequence composition, construct provenance and available training diversity. A benchmark must therefore evaluate unseen proteins, protect exposure history, and treat independent proteins rather than DNA observations as the statistical units.

## Results 1: Leakage-controlled benchmark

We evaluate canonical reverse-complement DNA classes under protein-cluster leave-one-out evaluation. All seven GSE237017 proteins are now development-exposed, so they are used for benchmark diagnosis and not independent validation. The corrected v0.6 M3 median Spearman is 0.1415; the sequence k-mer baseline is 0.2321. Per-protein points and seed stability are retained rather than hidden by a single median.

## Results 2: Sequence-conditioned transfer failure

Protein-conditioned models do not show a conditional advantage on the exposed cohort. Protein-shuffle prediction correlation is 1.0000 and target-shuffle correlation is 0.9969 for the corrected M3 diagnostic. This is conditioning collapse, not evidence that the model learned protein-specific specificity.

## Results 3: Representation repair and falsification

The v0.5 reverse-complement input averaging caused non-equivalent 7-mers to collide. v0.6 uses an injective RC-class representation and passes the synthetic conditionality benchmark. Ranking improves, but real-data conditioning remains collapsed. Dense supervision and local/residue attention pilots likewise fail to establish a conditional advantage. These are falsification results for specific implementation and scaling explanations, not proof that every architecture must fail.

## Results 4: Structure-method audit

DeepPBS and NA-MPNN are evaluated only where a legal protein-DNA structure exists. Both cover DBP35 and DBP48 (2/7): DeepPBS median Spearman 0.1593 and NA-MPNN -0.0407. Their PWM/PPM-derived scores are ranking proxies, not PBM E-scores or affinity. AF3/ContactSeek feasibility is documented but not promoted to a full landscape benchmark. These results do not justify a universal claim about structure-aware methods.

## Results 5: Failure modes and domain shift

The preregistered hard-case set contains 1,515 sequence-vs-experiment disagreements; 871 (57.5%) are unresolved by all current sequence models. Natural PBM and designed uPBM differ in k-mer length, score semantics, protein length/composition and construct provenance. These observations motivate protein-level diversity and domain/assay shift hypotheses without identifying a unique cause.

## Results 6: Preregistered independent validation

Reserved for preregistered independent validation. No external quantitative labels have been read. RFdiffusion3 is currently a metadata-only Tier 2 candidate; exact constructs, randomized-region schema and processed data availability remain unverified.

## Discussion

The strongest claim supported by the current evidence is limited to sequence-based approaches on the present exposed de novo DBP cohort. Large numbers of DNA observations do not substitute for independent protein diversity. The experiment-vs-experiment reference (approximately 0.5914 Spearman) is an empirical reproducibility context and must not be treated as a ceiling. Structure-aware generalization remains unresolved because coverage is incomplete.

## Limitations

The designed cohort contains seven proteins, all development-exposed. Structure methods cover only two. Natural-to-designed transfer uses unmatched assay and construct protocols. Hard-case categories depend on the frozen assay-derived labels. No independent confirmatory cohort is available yet.

## Methods and data availability

All frozen manifests, source hashes, scoring rules, exposure status and assay-tier endpoints are in `metadata/v0_9_external/`, `docs/v0_9_external/` and prior v0.5-v0.8 directories. The external runner refuses execution while `label_access` is false.
