# Decisions

1. Use RCSB PDB structural complexes as the first real paired DBP-DNA benchmark because they provide direct provenance and sequence data with minimal ambiguity.
2. Historical v0.1 rule, superseded in v0.2: select one protein entity and one DNA entity per PDB entry by choosing the longest protein and longest DNA chain when multiple candidates exist. This rule is retained only as an uncurated raw/interim collection heuristic.
3. Keep raw downloads under `data/raw/` and never overwrite them.
4. Use Python 3.13 in `.venv313` because Python 3.14 did not have fast wheel availability for the scientific stack.
5. Treat all current scores as `sequence-only proxy baseline` values.
6. Use chr22 for the genome-scan demo to keep the first pass small and reproducible.
7. Favor simple, inspectable retrieval logic over heavy external aligner setup for the demo stage.
8. Publish the project as a public GitHub repository under the existing `mcxmx` account using the stored Git credential, and push the current `qwqaq` history as the remote `main` branch.
9. For benchmark v0.2, treat RCSB PDB complexes as a structural cognate layer, not as quantitative specificity ground truth. A PDB structure alone sets `has_structural_cognate=True` but does not set `has_quantitative_specificity_ground_truth=True`.
10. For benchmark v0.2, keep the old longest-chain collection rule only as an uncurated raw/interim heuristic. Curated benchmark tables use biological mechanism annotation plus chain-contact evidence where available.
11. Use a small JASPAR CORE PFM-derived Layer C pilot for benchmark v0.2 because it provides accessible motif-backed quantitative DNA preference profiles with protein identity and UniProt provenance. These scores are stored as `jaspar_pfm_pwm_log2_odds_derived`, not as raw PBM or HT-SELEX enrichment.
12. For benchmark v0.3, use GEO GSE237017 processed uPBM 7-mer E-scores as the primary designed-DBP experimental specificity score. Store the consensus mean as `experimental_score_primary` and preserve raw/secondary score summaries separately.
13. Interpret v0.3 PBM E-scores only as per-protein experimental 7-mer specificity/enrichment scores. Do not treat them as Kd, binding free energy, binding probability, in vivo binding, or absolute cross-protein affinity.
14. The GSE237017 processed 7-mer files contain a primary 7-mer column and a reverse-complement companion 7-mer column. v0.3 explicitly expands both columns and verifies reverse-complement score consistency instead of silently merging or dropping one orientation.
15. Keep GSE237017 designed DBPs as an external/OOD designed-protein evaluation set for v0.4+ model development. Do not mix these seven designed proteins into natural DBP training splits unless the analysis is explicitly leave-one-designed-DBP-out.
16. For v0.3.1, treat 8192 source-data rows as reverse-complement equivalence classes and 16384 expanded sequences as oriented rows. Report 57344 protein-RC-class units as the primary independent sequence-unit count.
17. For v0.3.1, reproduce paper uPBM motif percentiles from Nature Source Data Extended Data Fig. 8 using 8192 RC-class rows, motif matching against either 7-mer column, rank-percentile E-scores within each replicate, and replicate averaging.
18. For v0.3.1, separate original design target, experimental assay target, and PBM evaluation motif. DBP48 is recorded as original design target I and assay/PBM reference target C.
19. For v0.3.1, define disagreement candidates as sequence-vs-experiment cases, not model failures, because no protein-conditioned model has been evaluated yet.
20. For future DNA-level splits, group by canonical reverse-complement equivalence class and never split a sequence and its reverse complement across train/test.
