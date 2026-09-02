# TODO

- Add additional processed experimental PBM/uPBM, raw array-level PBM where needed, HT-SELEX, and CIS-BP quantitative specificity sources beyond GSE237017 and the current JASPAR PFM-derived pilot.
- Decide whether each experimental score type can support ranking, regression, AUROC/AUPRC, or only descriptive analysis.
- Add protein-family and DNA-similarity split generation.
- Add leakage checks between structural cognate targets and experimental specificity data.
- Add calibration evaluation once true binding labels or calibrated quantitative assay scores are available.
- Train/evaluate `SimpleProteinConditionalBaseline` only after adding assay-matched natural PBM/uPBM training data; current v0.4 implementation is an untrained scaffold.
- Re-run DeepPBS in a supported Linux/container environment, or document it as not comparable if the official preprocessing cannot be reproduced.
- Expand structure-aware baseline coverage beyond the current DBP35/DBP48 NA-MPNN diagnostics.
- Keep GSE237017 designed DBPs as external/OOD test data unless explicitly running leave-one-designed-DBP-out experiments.
- Add assay-matched natural PBM/uPBM controls before making strong natural-to-designed OOD claims.
- Extend genome candidate retrieval from chr22 demo to a documented full-GRCh38 workflow.
