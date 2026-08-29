# TODO

- Add raw PBM / HT-SELEX / CIS-BP quantitative specificity sources beyond the current JASPAR PFM-derived pilot.
- Decide whether each experimental score type can support ranking, regression, AUROC/AUPRC, or only descriptive analysis.
- Add protein-family and DNA-similarity split generation.
- Add leakage checks between structural cognate targets and experimental specificity data.
- Add calibration evaluation once true binding labels or calibrated quantitative assay scores are available.
- Implement a real protein-conditioned scoring backend with `is_protein_conditioned=True`.
- Extend genome candidate retrieval from chr22 demo to a documented full-GRCh38 workflow.
