# V0.7 Scientific Assessment

## Four-state interpretation

The v0.5/v0.6 conditional models occupy **low-to-moderate performance + collapsed conditioning** on the exposed designed replay: v0.6 all-seven medians are M0 0.2829, M1 0.2016, M1c 0.1693, M2 0.1342 and M3 0.1415, while protein-shuffle median correlation is 1.0000 and target-shuffle median correlation is 0.9969. The sequence-only comparator is 0.2321. This prevents interpreting raw Spearman as proof of protein-specific specificity.

## Evidence and limits

There is sufficient evidence for a **benchmark-observed generalization gap under this strict unseen-protein protocol**: the gap is replicated across v0.5, v0.6, local/residue and dense-supervision diagnostics, and 871/1,515 hard cases remain unresolved by all current models. It is not yet sufficient to claim a universal biological impossibility. Seven independent designed proteins make a small-N artifact plausible; the simulation explicitly demonstrates that large DNA N cannot substitute for protein-level diversity, but it is not an empirical threshold.

## Formal claims now supported

1. Under the current GSE237017 design, strict unseen-protein ranking is modest and conditional heads do not show measurable protein sensitivity.
2. v0.5 RC averaging caused severe cross-class representation collisions; v0.6 repair improves developmental ranking but does not remove collapse.
3. Dense pair scaling and residue/local attention pilots did not establish a conditional advantage.
4. Natural PBM and designed uPBM differ in k-mer length, assay, score semantics and construct provenance; existing SimplePC transfer numbers are unmatched context.

## Hypotheses only

Natural-to-designed domain shift, assay/construct mismatch, insufficient protein diversity and inadequate sequence representations may each contribute. Their causal contributions cannot be separated with seven exposed proteins and unmatched natural constructs.

## Route decision

**Route B: benchmark / generalization-gap paper**, with a data-expansion prerequisite. Route A is not supported because no current conditional method passes real-data conditionality. Route C remains a requirement for a stronger causal claim: obtain a new independent designed-DBP cohort and preregister the frozen evaluation before reading labels. Do not present GSE237017 as Figure 5 external validation.

## Critical missing experiment

The single highest-value experiment is a genuinely independent, assay-matched designed-DBP cohort with substantially more than seven protein units, full 7-mer specificity landscapes, construct sequences and a preregistered frozen benchmark. This simultaneously tests whether the observed gap survives exposure control and whether protein-level sample size, rather than DNA-level pair count, is limiting.

## Strongest claim

**Under strict unseen-protein evaluation, existing sequence-based approaches show a reproducible benchmark-level generalization gap on the currently available de novo DBP cohort: representation repair improves ranking but does not restore protein- or target-dependent prediction, and the evidence is not yet sufficient for a new predictive-method claim or independent validation claim.**
