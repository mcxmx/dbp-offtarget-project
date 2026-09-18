# Submission Gap Analysis

This is a qualitative readiness assessment, not an acceptance-probability estimate. It is based on the frozen v1.3.4 evidence and assumes no new model training.

## Nucleic Acids Research

- **Novelty:** Strongest fit is the explicit separation of protein conditionality, position sensitivity, and nucleotide identity, with a prospective locked holdout and counterfactual analysis.
- **Methodological contribution:** Clear evaluation framework and reproducible forensic/provenance rules; variance allocation is useful if kept distinct from performance decomposition.
- **Dataset breadth:** Seven development proteins and three holdouts are informative but small.
- **Prospective validation:** Present and genuinely locked, but limited to three proteins.
- **External support:** SaMBA provides independent technical repeatability, while natural-TF structural prediction remains not evaluable.
- **Biological interpretation:** Appropriate for a methods/evaluation paper; claims about natural-TF generalization must remain limited.
- **Remaining weaknesses:** No competition replicate-level noise ceiling; conditionality collapse relies on cached historical diagnostics; no identity-aware model improvement is demonstrated.
- **Additional experiment:** A broader protein-level benchmark and replicate-resolved competition data would materially strengthen the submission, but are not prerequisites for a carefully scoped evaluation manuscript.

## Bioinformatics

- **Novelty:** Strong as a benchmark and diagnostic framework, especially the position-only counterfactual.
- **Methodological contribution:** Metrics, frozen split contract, variance allocation, and audit trail are concrete and reproducible.
- **Dataset breadth:** Likely the main concern; the manuscript should emphasize that the contribution is evaluation design rather than a universal model ranking.
- **Prospective validation:** A meaningful positive feature despite the small holdout.
- **External support:** SaMBA supports assay-independent technical repeatability but is not a matched benchmark.
- **Biological interpretation:** Keep conservative and method-specific.
- **Remaining weaknesses:** No new predictive baseline or identity-aware model is trained in this phase; readers may ask whether the counterfactual adds information beyond residual metrics.
- **Additional experiment:** Optional broader benchmark or a lightweight baseline trained on a defensible protein-level dataset would improve competitiveness, but should be a separate v1.4 project.

## PLOS Computational Biology

- **Novelty:** Good if framed around a general evaluation failure mode rather than a DeepPBS critique.
- **Methodological contribution:** The orthogonal axes and explicit counterfactual provide a conceptual contribution with reusable analysis code.
- **Dataset breadth:** Small sample size and source-specific designs limit broad biological claims.
- **Prospective validation:** Supports reproducibility of the main directional finding.
- **External support:** SaMBA is useful as a convergent technical control, with clear estimand limitations.
- **Biological interpretation:** Discussion should focus on what can and cannot be inferred from competition landscapes.
- **Remaining weaknesses:** Missing replicate-level competition data and unresolved natural-TF challenge.
- **Additional experiment:** Author-provided replicate records or an independent competition assay would be the highest-value addition.

## Nature Communications (stretch)

- **Novelty:** The evaluation concept is potentially broad, but current evidence is unlikely to support a high-impact generalization claim by itself.
- **Methodological contribution:** Reproducible decomposition and prospective design are strengths.
- **Dataset breadth:** The principal limitation is the small designed-DBP sample and lack of evaluable natural-TF prediction.
- **Prospective validation:** Positive but narrow.
- **External support:** SaMBA is orthogonal rather than matched, and does not close the biological validation gap.
- **Biological interpretation:** The paper would need a substantially broader protein-level benchmark, matched replicate data, or a demonstrated identity-aware model improvement.
- **Remaining weaknesses:** No causal intervention, no competition noise ceiling, and no broad natural-domain result.
- **Additional experiment:** A multi-protein, replicate-resolved, protein-held-out benchmark plus an identity-aware model would likely be required for this venue.

## Readiness decision

The evidence is sufficient to write and submit a scope-limited methods/evaluation paper to a computational or nucleic-acid-focused venue after normal manuscript polishing and author/source checks. It is not sufficient to claim universal natural-TF generalization or to justify a Nature Communications-level biological conclusion. The v1.4 identity-aware model gate remains closed until a defensible training set beyond the seven designed proteins and stronger competition-level identifiability evidence are available.
