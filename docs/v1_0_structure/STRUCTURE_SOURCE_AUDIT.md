# V1.0 Structure Source Audit

## Scope

This development-stage audit checked the local repository, the official `cjg263/dbp_design` checkout, Nature Supplementary Data 6, and RCSB PDB 8TAC. The primary source is [Glasscock et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12618268/), whose data-availability statement identifies Supplementary Data 6 as the PDB design-model archive and 8TAC as the DBP48 co-crystal. No v0.5-v0.9 file was modified. All seven GSE237017 proteins remain development-exposed.

## Finding

Supplementary Data 6 contains design-hit PDBs for DBP001, DBP003, DBP005, DBP006, DBP009, DBP035 and DBP048. In each recovered file, protein chain A matches the corresponding PBM protein sequence exactly and DNA chain B contains the intended target duplex/context. This changes the interpretation of the old 2/7 structure coverage: it was a resource-audit/pipeline omission, not evidence that five proteins lack structure or that structure methods are biologically inapplicable.

DBP48 also has experimental co-crystal [PDB 8TAC](https://www.rcsb.org/structure/8TAC). Its deposited construct includes terminal additions relative to the 65-aa PBM sequence, and the 21-nt crystal DNA assembly differs from the 14-bp PBM target. It is therefore a construct-sensitivity sanity check, not an exact assay-matched replacement.

The complete machine-readable chain, sequence, DNA, provenance and mismatch inventory is `metadata/v1_0_structure/structure_inventory.csv`. The downloaded archive hash and extraction provenance are recorded in `results/v1_0_structure/v1_0_run_manifest.json`.
