# V0.8 SOTA Results

## Structure diagnostics

DeepPBS: median Spearman 0.159271 across 2/7 proteins (DBP35 and DBP48 only). NA-MPNN is evaluated on the same two structure-covered proteins; its exact metrics are in `na_mpnn_metrics.csv`. Neither result is a seven-protein transfer estimate.

The experimental replicate reference is approximately 0.5914 Spearman (experiment-vs-experiment), a context reference rather than a mathematical ceiling. Sequence-only and v0.6 conditional rows remain exposed developmental comparators.

## Interpretation

The available structure rows do not establish that structure methods close the gap because coverage is only N=2 and both structures are development-exposed. AF3/ContactSeek feasibility does not create a full-landscape score. The strongest defensible claim remains limited to sequence-based approaches on the current cohort.

## Scientific judgment

1. The present evidence supports a reproducible sequence-model benchmark gap, not a universal computational impossibility.
2. DeepPBS and NA-MPNN do not currently provide evidence of a rescue; their partial medians are below the replicate reference and are underpowered.
3. The 871/1515 common failures cannot be declared structure-method failures because only a subset is structure-evaluable.
4. The paper title should say **sequence-based specificity prediction** unless a future, adequately covered structure benchmark changes the result.
5. No independent confirmatory cohort is currently qualified without label access and preregistration; the project is **BLOCKED ON INDEPENDENT CONFIRMATORY DATA**.
6. Route B remains the appropriate paper route, with data expansion as the key blocker.
