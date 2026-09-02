# v0.4.1 Summary Report

Date: 2026-09-02

## Natural PBM

- Final natural proteins: 57
- Coarse protein families: 7
- Species: 7
- Experiments/profile groups: 69
- Protein-DNA units: 1,875,072
- Assay: UniPROBE universal PBM, processed contiguous 8-mer E-score.
- K-mer length: 8 bp RC classes.
- Protein sequence completeness in final benchmark: 100%; construct sequence completeness remains 0%.
- Replicated protein/construct groups in downloaded metadata: 10; median Spearman 0.861 where available.
- 40% cluster split: train 39, validation 9, natural_test 9 proteins.

## Baselines

- Sequence-only designed median Spearman: 0.232.
- SimplePC natural held-out median Spearman: 0.301.
- SimplePC designed external median Spearman: 0.362.
- NA-MPNN diagnostic designed median Spearman: -0.041 over 2 proteins.
- DeepPBS: not run; Docker/WSL unavailable on this host.

## Failure Summary

SimplePC improves designed external median Spearman over the best sequence-only metric by 0.130, but remains 0.229 below the empirical replicate reference. It resolves 333/1,515 pre-registered sequence-vs-experiment disagreement candidates. DBP6 and DBP48 remain the clearest designed-protein failures for this baseline.
