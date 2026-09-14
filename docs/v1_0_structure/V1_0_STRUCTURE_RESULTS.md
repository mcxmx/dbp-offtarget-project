# V1.0 Structure-Conditioning Mechanism Pilot Results

## Status

This is a development/mechanism pilot on the seven GSE237017 DBPs. All seven proteins were development-exposed; none of these results is independent confirmation.

## Structure audit

The Supplementary Data 6 audit recovered design complex PDBs for all seven DBPs. Therefore v0.8's 2/7 structure coverage was a previous resource/pipeline audit gap.

## Frozen inference status

NA-MPNN specificity inference ran locally on all seven recovered design PDBs with the frozen `s_70114.pt` checkpoint. Its PPM-derived landscape metrics are in `results/v1_0_structure/nampnn_metrics.csv` and full landscapes are in `results/v1_0_structure/nampnn_results.csv`; median Spearman over evaluated proteins is 0.1247 (N=7). DeepPBS reproduction is recorded honestly in `deeppbs_results.csv`: the required torch_geometric and Linux 3DNA/X3DNA preprocessing stack is unavailable in this workspace, so the historical DBP35/DBP48 outputs are reused and DBP5/6/9 published applicability is not claimed as reproduced. This is an environment stop, not a biological failure.

### NA-MPNN per-protein ranking

| protein | status | Spearman |
|---|---|---|
| DBP1 | evaluated | -0.3110 |
| DBP3 | evaluated | -0.0774 |
| DBP5 | evaluated | 0.1839 |
| DBP6 | evaluated | 0.1247 |
| DBP9 | evaluated | 0.1043 |
| DBP35 | evaluated | 0.2420 |
| DBP48 | evaluated | 0.2675 |

## Conditionality

The NA-MPNN landscape correlation table is `protein_landscape_correlations.csv`. A structure model is considered protein-sensitive only descriptively from these between-protein landscapes; no arbitrary rescue threshold is imposed. The current output must not be interpreted as proof of transferability because the structures are design models and the cohort is exposed.

Across the 21 protein pairs, the median between-protein landscape Spearman correlation is 0.0059 (range -0.2941 to 0.7860 in the generated table). Thus the frozen structural input produces protein-dependent landscapes, while the median PBM ranking remains 0.1247.

## Contacts and controls

`interface_contact_features.csv` contains exploratory 4-A heavy-atom/base contact descriptors. `dbp48_design_vs_crystal.csv` compares the two DBP48 NA-MPNN landscapes descriptively; differences cannot be attributed to design-model error without matched construct and DNA. Rosetta and AF3/ContactSeek subsets are frozen but not run because the required runtimes are unavailable; their candidate manifests preserve the failure reasons.

### Design-interface descriptors

| protein | heavy contacts (4A) | base contacts (4A) | base-contacting residues |
|---|---:|---:|---:|
| DBP1 | 98 | 48 | 12 |
| DBP3 | 102 | 39 | 16 |
| DBP5 | 80 | 25 | 14 |
| DBP6 | 82 | 35 | 10 |
| DBP9 | 75 | 33 | 11 |
| DBP35 | 80 | 32 | 16 |
| DBP48 | 81 | 38 | 10 |

## H1-H5 judgment

H1/H2: explicit geometry is now demonstrably available as an input and NA-MPNN produces protein-specific PPMs, but a full DeepPBS reproduction and a fair cross-method conditionality comparison remain incomplete.

H3: no claim of reaching the experimental replicate reference is made from the partial/heterogeneous structure results.

H4: base-contact and register descriptors provide mechanistic hypotheses only (N=7).

H5: sequence-only failure cannot yet be attributed solely to missing geometry, because DeepPBS reproduction on DBP5/6/9 remains blocked by environment and no AF3/Rosetta diagnostic was executed.

**Scientific answer:** explicit interface geometry has not yet been shown to rescue the protein-conditioning failure. The correct v1.0 status is **INCONCLUSIVE / BLOCKED ON REPRODUCIBLE STRUCTURE-METHOD EXECUTION**, with the important correction that the former 2/7 coverage was an audit gap rather than biological unavailability.
