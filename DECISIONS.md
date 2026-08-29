# Decisions

1. Use RCSB PDB structural complexes as the first real paired DBP-DNA benchmark because they provide direct provenance and sequence data with minimal ambiguity.
2. Select one protein entity and one DNA entity per PDB entry, choosing the longest protein and longest DNA chain when multiple candidates exist.
3. Keep raw downloads under `data/raw/` and never overwrite them.
4. Use Python 3.13 in `.venv313` because Python 3.14 did not have fast wheel availability for the scientific stack.
5. Treat all current scores as `sequence-only proxy baseline` values.
6. Use chr22 for the genome-scan demo to keep the first pass small and reproducible.
7. Favor simple, inspectable retrieval logic over heavy external aligner setup for the demo stage.
8. Publish the project as a public GitHub repository under the existing `mcxmx` account using the stored Git credential, and push the current `qwqaq` history as the remote `main` branch.
9. For benchmark v0.2, treat RCSB PDB complexes as a structural cognate layer, not as quantitative specificity ground truth. A PDB structure alone sets `has_structural_cognate=True` but does not set `has_quantitative_specificity_ground_truth=True`.
10. For benchmark v0.2, keep the old longest-chain collection rule only as an uncurated raw/interim heuristic. Curated benchmark tables use biological mechanism annotation plus chain-contact evidence where available.
11. Use a small JASPAR CORE PFM-derived Layer C pilot for benchmark v0.2 because it provides accessible motif-backed quantitative DNA preference profiles with protein identity and UniProt provenance. These scores are stored as `jaspar_pfm_pwm_log2_odds_derived`, not as raw PBM or HT-SELEX enrichment.
