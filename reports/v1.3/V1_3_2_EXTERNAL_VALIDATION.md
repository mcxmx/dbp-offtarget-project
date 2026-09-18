# v1.3.2 External Validation

Status: **SaMBA technical-repeatability analysis complete; T227 row-level benchmark not recovered; v1.4 gate remains closed.**

## SaMBA technical repeatability

The primary canonical Watson-Crick analysis uses 10 source-validated sites across 7 TFs. It compares two deterministic spot-level technical estimates from the same official calibration arrays. The primary per-site results are in `results/v1_3/external_samba_canonical_decomposition.tsv`; TF-level summaries are in `results/v1_3/external_samba_decomposition_summary.tsv`.

| canonical metric | TF-level median | n TF |
|---|---:|---:|
| global Spearman | 0.947 | 7 |
| position-sensitivity Spearman | 0.881 | 7 |
| within-position residual Spearman | 0.823 | 7 |
| within-position pairwise accuracy | 0.814 | 7 |

For the pre-registered 100-seed secondary split stability analysis, the TF-level median identity residual Spearman had a 2.5%-97.5% split range of 0.781-0.874; pairwise accuracy had a corresponding range of 0.772-0.854. These are split-stability summaries over technical spot assignments, not biological confidence intervals.

The secondary mismatch analysis remains separate. Its TF-level medians were global 0.939, position 0.915, identity residual 0.882, and pairwise accuracy 0.803. These values are mechanistic/secondary evidence and are not pooled with canonical mutations.

## Interpretation boundary

The canonical technical split shows that identity-level effects are repeatable under the same SaMBA calibration assay. This makes a pure "identity is unobservable in any experiment" explanation less likely for this assay, but it does not establish reproducibility across biological experiments, assay contexts, or competition measurements. It also does not by itself prove that a DeepPBS identity failure is model-specific, because SaMBA and competition have different estimands and the designed competition workbook still lacks public replicate-level values.

The landscape variance fractions in the decomposition tables are descriptive allocations of experimental effect variance to position and residual components. They are not percentages of predictive correlation explained by position.

## T227 and v1.4 gate

The official T227 row-level table and published per-row predictions were not recovered. Therefore no natural-TF global/position/identity generalization claim is made. The external evidence currently supports strong SaMBA technical identity repeatability but lacks the natural-protein row-level method benchmark needed for the pre-registered v1.4 gate. No new GNN or other model was trained.

## Reproducible outputs

- Spot table: `data/processed/v1_3_2_samba_spots.parquet`
- Primary and secondary decomposition: `results/v1_3/external_samba_canonical_decomposition.tsv`, `external_samba_mismatch_decomposition.tsv`, `external_samba_decomposition_summary.tsv`
- Technical nulls: `results/v1_3/external_samba_permutation_nulls.tsv`, `external_samba_permutation_distributions.parquet`
- Figure source and render: `results/v1_3/figure3_external_identifiability.tsv`, `figures/v1_3/figure3_external_identifiability.png`, `figures/v1_3/figure3_external_identifiability.pdf`
