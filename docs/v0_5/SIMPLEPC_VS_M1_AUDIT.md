# SimplePC versus Current v0.5 M1 Audit

The prior v0.4.1 `SimpleProteinConditionalBaseline` and current v0.5 `M1`
numbers must not be interpreted as a head-to-head model comparison.

| Dimension | Prior SimplePC | Current v0.5 M1 |
| --- | --- | --- |
| Training proteins | Natural PBM benchmark | Designed DBP proteins in each LOCO training fold |
| Test regime | Designed external evaluation | Held-out designed protein cluster |
| Protein input | Frozen/compressed composition-style representation | Frozen ESM-2 t12 35M, 480-dimensional embedding |
| DNA input | Prior simple sequence feature representation | RC-symmetric one-hot 7-mer |
| Objective | Prior v0.4.1 regression/ridge-style protocol | Within-protein logistic pairwise ranking |
| Pair sampling | Prior protocol | 512 deterministic within-protein pairs |
| Evaluation unit | Prior designed external set | Protein first, 8,192 RC classes per protein |
| Target input | Not a target-conditioned model | M1 does not receive T |
| Leakage controls | v0.4.1 natural/designed protocol | Protein-cluster LOCO, no row-level split |

The most important difference is training regime: the prior result transfers
from a larger natural PBM collection, whereas current M1 learns its head from
only the other designed proteins in each fold. The comparison therefore mixes
training distribution, protein holdout, loss, representation, and sample
efficiency. It cannot identify a single causal reason for the numerical gap.

The current Phase 4 analysis deliberately does not run a bridge experiment.
Any bridge would be a new training diagnostic and would not alter the frozen
primary result.
