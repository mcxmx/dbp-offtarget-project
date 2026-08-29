# Decisions

1. Use RCSB PDB structural complexes as the first real paired DBP-DNA benchmark because they provide direct provenance and sequence data with minimal ambiguity.
2. Select one protein entity and one DNA entity per PDB entry, choosing the longest protein and longest DNA chain when multiple candidates exist.
3. Keep raw downloads under `data/raw/` and never overwrite them.
4. Use Python 3.13 in `.venv313` because Python 3.14 did not have fast wheel availability for the scientific stack.
5. Treat all current scores as `sequence-only proxy baseline` values.
6. Use chr22 for the genome-scan demo to keep the first pass small and reproducible.
7. Favor simple, inspectable retrieval logic over heavy external aligner setup for the demo stage.

