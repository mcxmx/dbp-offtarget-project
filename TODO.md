# TODO

- Add additional processed experimental PBM/uPBM, raw array-level PBM where needed, HT-SELEX, and CIS-BP quantitative specificity sources beyond GSE237017 and the current JASPAR PFM-derived pilot.
- Decide whether each experimental score type can support ranking, regression, AUROC/AUPRC, or only descriptive analysis.
- Add protein-family and DNA-similarity split generation.
- Add leakage checks between structural cognate targets and experimental specificity data.
- Add calibration evaluation once true binding labels or calibrated quantitative assay scores are available.
- Re-run DeepPBS in a supported Linux/Docker/WSL environment using `external/deeppbs/run_official_example_linux.sh`, then evaluate structurally available designed DBPs.
- Replace v0.4.1 amino-acid 3-mer Jaccard split proxy with MMseqs2 or CD-HIT clustering.
- Curate UniPROBE experimental construct sequences; current v0.4.1 natural benchmark uses full-length UniProt references marked `sequence_match_to_assay=false`.
- Expand natural PBM benchmark beyond the current 57 train-ready proteins after sequence/provenance curation.
- Expand structure-aware baseline coverage beyond the current DBP35/DBP48 NA-MPNN diagnostics.
- Keep GSE237017 designed DBPs as external/OOD test data unless explicitly running leave-one-designed-DBP-out experiments.
- Add assay/k-mer-length matched natural uPBM controls before making strong natural-to-designed OOD claims.
- Extend genome candidate retrieval from chr22 demo to a documented full-GRCh38 workflow.
