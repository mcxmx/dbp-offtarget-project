# v0.4.1 K-mer Length Harmonization

## Natural PBM

The v0.4.1 natural PBM source uses UniPROBE contiguous 8-mer E-scores. Each source row contains an 8-mer and its reverse complement plus one E-score, so the natural sequence unit is an 8-mer reverse-complement equivalence class.

## Designed uPBM

The v0.3/v0.3.1 designed DBP benchmark from GSE237017 uses processed 7-mer PBM E-scores. Its independent DNA units are 7-mer reverse-complement classes.

## Score Compatibility

Both sources are PBM-derived enrichment/specificity scores, but their exact processing pipelines are not guaranteed identical. Scores are used for per-protein ranking, not cross-protein absolute affinity.

## Modeling Implication

Future protein-conditioned baselines must support variable-length DNA inputs or use a predeclared length harmonization strategy. v0.4.1 does not crop natural 8-mers into 7-mers and does not pad designed 7-mers into 8-mers.

## Confounding Risk

Natural held-out versus designed external performance can be confounded by k-mer length, array design, score processing, protein construct differences, and assay protocol. A drop on designed DBPs cannot by itself be interpreted as pure biological OOD.
