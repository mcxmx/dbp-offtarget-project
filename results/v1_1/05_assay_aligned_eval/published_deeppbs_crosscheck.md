# Published DeepPBS cross-check

The paper reports qualitative agreement between predicted DeepPBS PWMs and competition specificity for DBP5, DBP6, DBP9 and DBP35. The v1.1 quantitative analysis uses independently fixed assay positions and `delta_logP` without optimizing register, orientation, epsilon, sign, or any experimental metric.

| published protein | mutations | v1.1 Spearman |
|---|---:|---:|
| DBP005 | 42 | 0.3492 |
| DBP006 | 42 | 0.7448 |
| DBP009 | 42 | 0.5655 |
| DBP035 | 42 | 0.3584 |

These values are a quantitative diagnostic cross-check, not a parameter-selection target. The source XLS provides the two-replicate heatmap means, not separate replicate columns. Any disagreement with the published qualitative display remains interpretable in light of register, strand, construct, and assay-direction checks documented in the source and mapping reports.
