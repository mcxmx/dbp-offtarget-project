# NA-MPNN training-overlap audit

Exact designed-protein training sequences are not distributed in the checked-out repository; homolog-level overlap is therefore UNDETERMINED.

## Explicit hits
- `external/nampnn/NA-MPNN/splits/design_valid.json` contains `8tac`: POSSIBLE/CONFIRMED overlap requires interpretation of split semantics.
- `external/nampnn/NA-MPNN/splits/specificity_valid.json` contains `8tac`: POSSIBLE/CONFIRMED overlap requires interpretation of split semantics.

DBP48/8TAC was historically flagged in local NA-MPNN split files; the v1.1 design-model runs use Supplementary Data 6 design PDBs, not 8TAC. This does not establish a clean homolog-level separation.
