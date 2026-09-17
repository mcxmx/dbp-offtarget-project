# v1.3.1 DBP35 to DBP35opt Analysis

## DATA-AVAILABILITY AMENDMENT

The frozen v1.3 evaluation contract pre-registered the real-perturbation endpoint but marked it unavailable. Official Extended Data Fig. 9 source data now make the experimental landscape component evaluable. This is a data-availability amendment, not a metric amendment: effect direction, `mean_abs_effect` position sensitivity, and within-position centering remain unchanged.

## Recovered data

DBP35opt is DBP35 K18V/R33N/P42Q. The official panel 9e sheet contains 42 matched substitutions across 14 positions. It provides a published summary column but no replicate-valued columns. DBP35opt used 20 nM biotinylated target and 160 nM competitor, whereas the standard DBP35 competition landscape used 1 uM target and 8 uM competitor. The comparisons below are therefore continuous descriptive comparisons with an assay-context confound, not an isolated causal protein-mutation effect.

## Experimental decomposition

- Global landscape conservation: Spearman 0.736; Pearson 0.847.
- Position-sensitivity conservation: Spearman 0.776.
- Within-position identity conservation: residual Spearman 0.184.
- Mean absolute mutation-landscape shift: 0.357 normalized PE/FITC units.
- Mean absolute identity-residual shift: 0.106.
- Positive continuous sensitivity change at positions: 1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14.
- Negative continuous sensitivity change at positions: 2, 6.

No threshold was chosen and no binary gained/lost-specificity claim is made. Exact mutation and position values are in `results/v1_3/dbp35opt_mutation_shift.tsv` and `results/v1_3/dbp35opt_position_shift.tsv`.

## Model endpoint

The predicted landscape-shift endpoint remains `NOT_EVALUABLE_NO_DBP35OPT_MODEL_PREDICTION`: the recovered workbook supplies experiment, not a DBP35opt structure/model prediction. It is therefore not yet possible to say whether DeepPBS detects this real protein perturbation.
