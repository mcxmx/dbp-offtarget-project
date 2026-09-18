# v1.3.3 SAMPDI Leakage Audit

Status: **SECONDARY / NOT USED FOR PRIMARY CHALLENGE**.

The v1.3.2 T227 audit remains frozen: the official SAMPDI-3D paper and
supplement expose aggregate benchmark values but no auditable row-level T227
mutation table or per-row SAMPDI/FoldX predictions. The v1.3.3 natural-TF
challenge therefore does not reopen T227 recovery and does not use
SAMPDI-3D/3Dv2 predictions on SaMBA.

Because no row-level source or training-overlap manifest was recovered, exact
TF/site/mutation overlap with SaMBA cannot be excluded. Any future SAMPDI
result on SaMBA must be labelled `DESCRIPTIVE / NON-INDEPENDENT` unless the
training data and exact overlap audit establish independence. It cannot replace
the frozen official DeepPBS primary challenge.
