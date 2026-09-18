# v1.3.1 Competition Replicate Analysis

## Source audit result

The official publisher workbook for Extended Data Fig. 3 contains one `Median PE/FITC (Normalized)` column per protein and no replicate identifier or replicate-valued columns. The downloadable workbook therefore does not contain the individual replicate experiments mentioned in the article caption. `data/processed/v1_3_competition_mutations_replicates.parquet` is an intentionally empty schema-only table; no replicate values were inferred from the published mean or figure pixels.

## Published-mean parity

The newly downloaded publisher workbook was joined to the repository's prior copy by protein, position, wild-type base, and mutant base. Across 414 mutations, maximum absolute error was 0 and median absolute error was 0. This verifies source-file parity, not replicate reconstruction.

## Noise ceiling status

Mutation-level, position-sensitivity, and within-position identity replicate agreement are `NOT_EVALUABLE_SOURCE_DATA_MEAN_ONLY`. The requested output `results/v1_3/replicate_noise_ceiling.tsv` retains all ten proteins with explicit missing metrics. Consequently, current data cannot determine whether weak identity-level model performance reflects a model-specific failure or an identity-level experimental reproducibility ceiling.

## PBM interpretation

`results/v1_3/assay_reproducibility_matrix.tsv` places PBM-versus-competition values beside explicit unavailable competition R1-versus-R2 cells at GLOBAL, POSITION, and IDENTITY levels. PBM disagreement remains a mixture of assay estimand, biological/context, and technical differences; technical variance is not separately identifiable without the missing competition replicates.

This evaluability boundary is shown without imputation in `figures/v1_3/Figure3_reproducibility_ceiling.png` and PDF, using `results/v1_3/figure3_reproducibility_data.tsv` as source data.
