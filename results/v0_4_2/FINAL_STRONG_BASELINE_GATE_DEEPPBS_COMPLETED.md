# DeepPBS-Completed Baseline Gate

The prior v0.4.2 `WAIT` status was caused by an unevaluated DeepPBS runtime,
not by a model result. This completion artifact supersedes that interpretation
without deleting the historical report.

## Gate

**WAIT FOR STRONGER STRUCTURE-COVERED BASELINE**

DeepPBS is now technically integrated and has a valid two-protein diagnostic
result, but the structure coverage is only **2/7**. The available
result does not justify calling DeepPBS a complete seven-protein strong
baseline or claiming designed-DBP generalization. The project may use it as a
structure-aware diagnostic comparator while pursuing additional public
designed complexes or a separately labeled predicted-structure sensitivity
study.

## Scientific conclusion

DeepPBS can be evaluated against the experimental PBM landscape when a
protein-matched protein-DNA structure is available. This run does not show
that structure-aware scoring solves the designed specificity problem: the
experimental reference is `0.591` median Spearman, while the current
DeepPBS coverage is only two proteins. The remaining five missing structures
are a coverage limitation, not missing predictions filled with zero.

Do not call DBP35/DBP48 a clean zero-shot result: exact overlap was not found in
the checked manifests, but homolog-level overlap remains unresolved.
