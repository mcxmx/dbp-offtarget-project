# Proposed Model Requirements V2

Date: 2026-09-02

These are conditional requirements, not an implemented architecture.

## Requirements Driven by v0.4.1

1. Support variable-length DNA k-mers at minimum 7-mer and 8-mer inputs without cropping/padding that changes the biological task.
2. Use protein information beyond amino-acid composition, because the SimplePC composition ridge still leaves a large gap to replicate agreement.
3. Preserve per-protein ranking evaluation and RC canonicalization as first-class constraints.
4. Explicitly model natural-to-designed generalization without conflating assay shift with biological OOD.
5. Include structure-aware baseline comparison once DeepPBS is runnable in Linux/Docker.

The final model should not be started until DeepPBS is fairly reproduced or formally ruled out as non-evaluable for the designed benchmark.
