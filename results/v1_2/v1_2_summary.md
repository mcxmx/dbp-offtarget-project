# v1.2 summary

This is a frozen-input robustness and assay-ceiling audit. It reads v1.1 artifacts without changing them; all seven GSE237017 proteins remain development-exposed.

## Primary results
- v1.1 reproduction: all seven DeepPBS local correlations passed at 1e-12 tolerance; primary-six median = 0.3538.
- Position-preserving null: 1/6 primary proteins exceed their null 97.5th percentile; see `02_permutation_nulls/permutation_summary.tsv`.
- DeepPBS primary-six local median = 0.3538; PBM experimental local median = 0.3636. These headline values use different mutation coverage; the matched 21-mutation comparison is in the assay triangle.
- Position-level DeepPBS Spearman values: DBP001 0.4615, DBP003 0.2582, DBP005 0.3846, DBP006 0.8022, DBP009 0.7495, DBP035 0.1516, DBP048 -0.0462
- Within-position pairwise ordering accuracy: primary-six median = 0.5064 (chance reference 0.5); top-deleterious-base exact agreement median = 0.2967 (chance reference 1/3).
- DeepPBS large-effect recovery: primary-six median AP = 0.5896.
- Leave-one-position-out primary ranges: 0.0574 to 0.7725; individual summaries are in `10_position_robustness/summary.tsv`.

## Per-protein DeepPBS local correlations

| protein | rho | position rho | within-position pairwise accuracy | top-base exact | mapping |
|---|---:|---:|---:|---:|---|
| DBP001 | 0.2271 | 0.4615 | 0.3846 | 0.3077 | PRIMARY |
| DBP003 | 0.1824 | 0.2582 | 0.5128 | 0.4615 | PRIMARY |
| DBP005 | 0.3492 | 0.3846 | 0.5000 | 0.2857 | PRIMARY |
| DBP006 | 0.7448 | 0.8022 | 0.4286 | 0.2857 | PRIMARY |
| DBP009 | 0.5655 | 0.7495 | 0.5714 | 0.2857 | PRIMARY |
| DBP035 | 0.3584 | 0.1516 | 0.6667 | 0.5714 | PRIMARY |
| DBP048 | 0.0498 | -0.0462 | 0.4524 | 0.3571 | AMBIGUOUS / SENSITIVITY_ONLY |

## Baseline and assay interpretation
- The position-only PPM-information baseline tests position sensitivity without mutant-base identity. The fixed transition indicator is a non-fitted, non-mechanistic mutation-type control. The uniform baseline is correctly undefined.
- PBM-vs-competition primary-six median = 0.3636; values range from -0.4221 to 0.8133.
- Exploratory across-protein association of assay agreement and DeepPBS performance: Spearman = 0.0286 (N=6); this is not an assay ceiling ratio.
- DeepPBS published four-protein descriptive subset median = 0.4620; this reuses the same exposed proteins and assay.
- NA-MPNN comparison is limited to the explicitly named frozen global-landscape-derived local score; it is not a native NA-MPNN PPM.
- Directionality is not applicable because the competition XLS does not define WT-relative signed improve/reduce effects.

## R1-R5 interpretation
- R1 genuine local structure signal: NOT ESTABLISHED; permutation excess and within-position results must be read protein by protein, without a population claim.
- R2 mostly position-sensitivity signal: SUPPORTED AS THE DOMINANT OBSERVED COMPONENT; only DBP035 exceeded the position-preserving null, while primary-six within-position pairwise accuracy was near its 0.5 chance reference.
- R3 mostly assay-shared signal: INCONCLUSIVE; PBM-vs-competition agreement is heterogeneous and no invalid normalization by assay rho is used.
- R4 heterogeneous/protein-specific success: SUPPORTED DESCRIPTIVELY; DBP006 and DBP009 drive much of the positive DeepPBS local median, while DBP001/003 are weaker.
- R5 fragile local result: NOT INDICATED BY SIGN CROSSING IN PRIMARY LOO; mutation jackknife and per-position ranges quantify influence without inferential claims.

## Boundary
The result does not establish transferable specificity, does not make the seven proteins independent validation, and does not represent all structure-based methods. The most important unresolved limitation remains an independent assay-matched designed-DBP cohort with frozen analysis before label access.
