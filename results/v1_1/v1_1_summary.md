# v1.1 summary

This is a development-stage assay-alignment audit on the seven development-exposed GSE237017 DBPs. No v0.x or v1.0 output was overwritten.

## Findings
1. DeepPBS design-model coverage/QC: 7/7 frozen NPZ predictions validated.
2. DeepPBS local competition performance: primary six-protein median Spearman 0.3538; DBP048 is sensitivity-only because sequence C differs from the design sequence.
3. NA-MPNN local performance: primary six-protein median Spearman 0.1487; every value is explicitly `NA-MPNN frozen global landscape-derived local mutation score`, not a native PPM result.
4. PBM-vs-competition agreement: primary six-protein median Spearman 0.3636; this compares experimental assays, not a model.
5. Global PBM comparison is retained in `local_global_comparison.tsv`; individual proteins, not pooled sequence counts, are the statistical units.
6. Experimental additive-PWM capacity at Hamming distance >=2 remains median 0.1294.
7. Local-to-global distance strata remain descriptive and are not used to select a register.
8. Global conditionality remains protein-specific: DeepPBS median pairwise correlation 0.2364; NA-MPNN 0.0059.
9. Training overlap remains UNDETERMINED because full checkpoint training sequences are unavailable.

## Per-protein local/global table

| protein | DeepPBS local | DeepPBS global PBM | NA-MPNN derived local | NA-MPNN global PBM | PBM vs competition | mutations | mapping |
|---|---:|---:|---:|---:|---:|---:|---|
| DBP001 | 0.2271 | -0.0661 | 0.2286 | -0.3110 | -0.4221 | 39 | PRIMARY |
| DBP003 | 0.1824 | 0.0083 | 0.3494 | -0.0774 | 0.1688 | 39 | PRIMARY |
| DBP005 | 0.3492 | -0.0017 | 0.3104 | 0.1839 | 0.1481 | 42 | PRIMARY |
| DBP006 | 0.7448 | 0.0357 | 0.0675 | 0.1247 | 0.8133 | 42 | PRIMARY |
| DBP009 | 0.5655 | 0.0263 | 0.0688 | 0.1043 | 0.6831 | 42 | PRIMARY |
| DBP035 | 0.3584 | 0.0249 | 0.0052 | 0.2420 | 0.5584 | 42 | PRIMARY |
| DBP048 | 0.0498 | 0.0800 | -0.0455 | 0.2675 | 0.4818 | 42 | AMBIGUOUS |

## A/B/C/D classification
A: supported descriptively. DeepPBS local median (0.3538) is above its global PBM median (0.0166); the NA-MPNN-derived local median (0.1487) is only modestly above its global median (0.1145). This is not an independent or population-level result.
B: not supported as a blanket claim. Local performance is heterogeneous: DeepPBS spans 0.1824 to 0.7448 among primary proteins, while the derived NA-MPNN scores are mostly weak.
C: partially supported as an interpretation constraint, not as a universal assay-disagreement claim. PBM-vs-competition primary rho spans -0.4221 to 0.8133 (median 0.3636), so assay agreement is protein-dependent.
D: supported as a representation-capacity diagnostic because the experimental additive model is strong locally but degrades outside Hamming-1.

Overall status: DEVELOPMENT-EXPOSED / ASSAY-ALIGNED DIAGNOSTIC COMPLETE. This does not constitute independent validation or a universal structure-method conclusion.
