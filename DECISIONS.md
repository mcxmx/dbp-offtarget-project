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
21. For v0.4, freeze v0.3.1 as the fixed experimental benchmark and write all new outputs under `v0_4` paths. Do not modify v0.3.1 benchmark files.
22. For v0.4, evaluate baselines as per-protein DNA ranking tasks on canonical RC classes, then macro-summarize. Do not pool all 57,344 rows into a primary correlation.
23. For v0.4, treat DeepPBS as not fairly evaluable in the current Windows environment unless the official preprocessing stack can be run reproducibly. Missing DeepPBS predictions are stored as not evaluable, never as zero scores.
24. For v0.4, map NA-MPNN specificity-mode PPM outputs to 7-mer rankings only as a documented diagnostic task-harmonization layer: best contiguous structure-window log probability over both RC orientations.
25. For v0.4, classify DBP48/8TAC NA-MPNN results as diagnostic rather than zero-shot because `8tac` appears in NA-MPNN split files.
26. For v0.4, keep `SimpleProteinConditionalBaseline` as an untrained Tier 1 scaffold until an assay-matched natural PBM/uPBM training source is available. Do not train it on random rows from the seven designed DBPs as a main generalization result.
27. For v0.4, use empirical replicate agreement only as an assay reproducibility reference, not as a strict theoretical upper bound.
28. For v0.4, set the new-model gate to `CONDITIONAL GO`: the sequence-only gap is clear, but structure-aware strong baseline coverage remains incomplete.
29. For v0.4.1, use UniPROBE processed contiguous 8-mer E-scores as the primary natural PBM training source. Keep natural 8-mer PBM and designed 7-mer uPBM as separate assay/task layers.
30. For v0.4.1, include only complete 8-mer reverse-complement class profiles with conservative protein sequence recovery in the train-ready natural benchmark. Incomplete profiles and unclear/fusion labels are retained in QC metadata but excluded from training.
31. For v0.4.1, use full-length UniProt reference sequences only as conservative protein features and mark `sequence_match_to_assay=false` until assay construct sequences are curated.
32. For v0.4.1, use a fast amino-acid 3-mer Jaccard proxy only for split hygiene. It is not a formal homology analysis and should be replaced by MMseqs2/CD-HIT in a Linux runtime.
33. For v0.4.1, train `SimpleProteinConditionalBaseline_composition_ridge` only on natural PBM train proteins and select hyperparameters only on natural validation proteins. Designed DBPs remain an external test set.
34. For v0.4.1, do not claim DeepPBS performance because Docker and WSL are unavailable on this host. Record official provenance and runnable Linux wrapper instead of substituting a reimplementation.
35. For v0.4.1, set the final gate to `WAIT FOR STRONGER BASELINE`: SimplePC improves over sequence-only, but DeepPBS still needs a fair supported-runtime reproduction before proposed model development.
36. For v0.5, define `primary_target` as the independently reported original design target. Never infer it from PBM top motifs, PBM E-scores, or test-set optima; preserve experimental assay/reference DNA and PBM-derived motif as separate metadata fields.
37. For v0.5, use the existing four protein-sequence clusters for the primary leave-one-cluster-out split. Use connected components over protein-cluster, original-target, assay-target, and motif links as a strict sensitivity split.
38. For v0.5, treat canonical reverse-complement equivalence classes as the DNA unit and prohibit random protein-7-mer row splits.
39. For v0.5, a future target-conditioned model must be non-separable in `(P,T,D)`. Subtracting a candidate-independent `S(P,T)` constant from `S(P,D)` cannot change the within-protein ranking.
40. For v0.5 Phase 2 smoke training, use deterministic within-protein logistic ranking pairs with seed 42 and fixed 40% easy, 35% medium, and 25% hard rank-difference strata. The sampler was corrected to match this declared protocol before the final smoke run; this changes only the smoke-training pair composition, not benchmark data or split definitions.
41. For v0.5 Phase 3, classify `protein_cluster_loco_fold_1` (DBP1/DBP3) as development-exposed because it was used for engineering smoke training. Keep folds 2-4 as the previously unseen primary subset, while still reporting all seven proteins.
42. For v0.5 Phase 3, freeze evaluation seeds to 17, 29, and 43 before running unseen folds. Aggregate seed results within each protein/model before calculating macro medians; seed variation is not biological uncertainty.
43. For v0.5 Phase 3, retain the fixed M0/M1/M1c/M2/M3 architectures, target definitions, split manifests, and pair sampler. The complete primary result is evidence for later hard-case analysis, not a final model GO/NO-GO decision.
44. For v0.5 paired model comparisons, use the per-protein seed-mean deltas and their macro median as the primary contrast. Do not infer paired improvement from the difference between unpaired macro medians, because the median of paired differences is not generally equal to the difference of medians.
