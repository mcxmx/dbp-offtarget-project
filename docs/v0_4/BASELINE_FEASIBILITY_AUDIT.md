# v0.4 Baseline Feasibility Audit

Audit date: 2026-09-02

## Benchmark Held Fixed

Input benchmark is v0.3.1: 7 designed DBPs, 57,344 protein-RC-class experimental units, primary score `experimental_escore_consensus` from processed GSE237017 uPBM E-scores. v0.3.1 files are not modified.

## DeepPBS

- Official code: https://github.com/timkartar/DeepPBS
- Local path: `external/deeppbs/DeepPBS`
- Commit: `8bfb211dd67f02877841f6f33aa493ddf7daedf9`
- Paper DOI from README: 10.1038/s41592-024-02372-w
- Official output: predicted DNA position weight matrix / PWM from a protein-DNA complex structure.
- Required input: biological protein-DNA assembly as PDB/CIF, then DeepPBS preprocessing into graph/shape features.
- Structure requirement: yes. The method is not a raw `(protein sequence, DNA 7-mer)` scorer.
- Fixed cognate structure to alternative DNA ranking: only possible after a justified PWM-to-7mer mapping; it does not directly predict PBM E-score.
- Designed DBP overlap: no direct GSE237017/DBP ID or 8TAC hit in checked DeepPBS fold files; full homolog overlap cannot be excluded without the complete training-set sequence list.
- v0.4 execution status: not evaluated. The official preprocessing chain is Linux-oriented and requires 3DNA/Curves-style processing plus `torch_geometric`, `freesasa`, and related dependencies. Current Windows environment is not a fair execution target.

## NA-MPNN

- Official code: https://github.com/baker-laboratory/NA-MPNN
- Local path: `external/nampnn/NA-MPNN`
- Commit: `9fabc2482092b725e067969fba21297a806b6fda`
- Official output: predicted PPM over residue types for nucleic-acid positions from a protein-nucleic-acid structure.
- Required input: protein-DNA/RNA structure.
- Structure requirement: yes. The specificity mode is not a structure-free `(protein sequence, DNA 7-mer)` scorer.
- Fixed cognate structure to alternative DNA ranking: possible only by mapping PPM log-probabilities over a structure-defined DNA window into candidate 7-mer scores. This is a derived task-harmonization layer.
- Designed DBP overlap: `8tac` appears in NA-MPNN `design_valid.json` and `specificity_valid.json`; DBP48/8TAC is therefore not zero-shot. No DBP35/GSE237017 direct hit was found in checked split files, but homolog-level risk is unresolved.
- v0.4 execution status: official inference ran successfully on the bundled 1am9 example, on `external/dbp_design/2b_design_mpnn/DBP035.pdb`, and on RCSB `8TAC`.

## TransBind

TransBind is not executed in v0.4. It remains a task-comparison item until its input/output can be fairly mapped to this processed uPBM 7-mer ranking benchmark.

## Structure Availability Summary

| protein_id | structure_type | experimental_or_predicted | pdb_id | local_file | source_url | paper_reference | bound_dna_sequence | bound_dna_length | protein_chain | dna_chain | protein_length_in_structure | structure_confidence | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DBP1 | not found in current official checkout | unknown/not available | nan | nan | https://github.com/cjg263/dbp_design | Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4 | nan | nan | nan | nan | nan | none | No protein-DNA complex structure/model for this designed DBP was found in the checked official dbp_design repository. Structure-aware baselines are not evaluated for this protein in v0.4. |
| DBP3 | not found in current official checkout | unknown/not available | nan | nan | https://github.com/cjg263/dbp_design | Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4 | nan | nan | nan | nan | nan | none | No protein-DNA complex structure/model for this designed DBP was found in the checked official dbp_design repository. Structure-aware baselines are not evaluated for this protein in v0.4. |
| DBP5 | not found in current official checkout | unknown/not available | nan | nan | https://github.com/cjg263/dbp_design | Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4 | nan | nan | nan | nan | nan | none | No protein-DNA complex structure/model for this designed DBP was found in the checked official dbp_design repository. Structure-aware baselines are not evaluated for this protein in v0.4. |
| DBP6 | not found in current official checkout | unknown/not available | nan | nan | https://github.com/cjg263/dbp_design | Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4 | nan | nan | nan | nan | nan | none | No protein-DNA complex structure/model for this designed DBP was found in the checked official dbp_design repository. Structure-aware baselines are not evaluated for this protein in v0.4. |
| DBP9 | not found in current official checkout | unknown/not available | nan | nan | https://github.com/cjg263/dbp_design | Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4 | nan | nan | nan | nan | nan | none | No protein-DNA complex structure/model for this designed DBP was found in the checked official dbp_design repository. Structure-aware baselines are not evaluated for this protein in v0.4. |
| DBP35 | designed complex model | predicted/theoretical Rosetta model | nan | external/dbp_design/2b_design_mpnn/DBP035.pdb | https://github.com/cjg263/dbp_design | Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4 | GCAGATCTGCACATCGATGTGCAGATCTGC | 30 | A | B | 63 | medium | Only readily available designed-DBP complex model found in the official dbp_design repository checkout. Header marks EXPDTA THEORETICAL MODEL / Rosetta. DNA chain B contains target-B duplex context; NA-MPNN featurizer masks 29 of 30 DNA residues. |
| DBP48 | experimental protein-DNA complex | experimental X-ray diffraction | 8TAC | data/raw/rcsb/mmcif/8TAC.cif | https://www.rcsb.org/structure/8TAC | Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4; PMID 40940539 | ACCTGACGCGA;TTCGCGTCAGG | 21 | B | C;D | 66 | high | RCSB mmCIF reports X-RAY DIFFRACTION at 2.34 A and chains A/B as the designed protein entity; v0.2 curation selected protein chain B and DNA chain C. NA-MPNN is run on the full biological assembly, so this is a structure-level diagnostic, not isolated single-chain scoring. |

## Main Feasibility Conclusion

v0.4 can fairly evaluate sequence-only baselines on all seven proteins and partial NA-MPNN structural-PPM baselines on DBP35 and DBP48. DeepPBS is recorded as not evaluable in this environment. DBP48/8TAC has explicit NA-MPNN split overlap risk and is not zero-shot.
