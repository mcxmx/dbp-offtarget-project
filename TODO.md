# TODO

- Add additional processed experimental PBM/uPBM, raw array-level PBM where needed, HT-SELEX, and CIS-BP quantitative specificity sources beyond GSE237017 and the current JASPAR PFM-derived pilot.
- Decide whether each experimental score type can support ranking, regression, AUROC/AUPRC, or only descriptive analysis.
- Add protein-family and DNA-similarity split generation.
- Add leakage checks between structural cognate targets and experimental specificity data.
- Add calibration evaluation once true binding labels or calibrated quantitative assay scores are available.
- Implement a v0.4 protein-conditioned baseline backend with `is_protein_conditioned=True`, starting with lightweight/non-neural ranking models before deep learning.
- Keep GSE237017 designed DBPs as external/OOD test data unless explicitly running leave-one-designed-DBP-out experiments.
- Add assay-matched natural PBM/uPBM controls before making strong natural-to-designed OOD claims.
- Extend genome candidate retrieval from chr22 demo to a documented full-GRCh38 workflow.
