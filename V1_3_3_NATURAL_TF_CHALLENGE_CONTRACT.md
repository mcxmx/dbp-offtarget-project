# v1.3.3 Natural-TF Structural Model Challenge Contract

Status: **FROZEN BEFORE NATURAL-TF MODEL/LABEL PERFORMANCE INSPECTION**, 2026-09-18.

This contract extends the frozen v1.3/v1.3.2 decomposition without changing any
primary endpoint. It defines a prospective challenge of the official DeepPBS
ensemble against the canonical Watson-Crick subset of the Afek et al. SaMBA
technical-repeatability dataset.

## 1. Eligible labels

- Input labels are `data/processed/v1_3_2_samba_spots.parquet`.
- The primary endpoint includes only rows with
  `perturbation_type == canonical_watson_crick`.
- Mismatch rows are excluded from the primary challenge and may appear only in a
  separately labelled secondary analysis.
- One experimental label per `(tf_id, binding_site_id, position, wt_basepair,
  perturbation)` is constructed from the already frozen source-defined
  `published_fold_change`: `log2(published_fold_change)`. No new normalization,
  centering, median/mean choice, or label-dependent transformation is allowed.
- The primary label set is additionally restricted to the source-integrity-valid
  canonical sites carried forward by v1.3.2:
  `Cbf1`, `Egr1`, `Ets1_site_1`, `Ets1_site_3`, `GR`, `Max_site_1`,
  `Max_site_2`, `TBP_site_1`, `TBP_site_2`, and `p53`. Ets1_site_2 and
  Ets1_site_4 remain auditable raw sites but are excluded from primary labels
  because their published medians failed the frozen source reconstruction.
- Larger values mean stronger relative binding. The DeepPBS score keeps the
  previously frozen sign: `log(P_mutant + 1e-12) - log(P_wild_type + 1e-12)`.

## 2. Structure eligibility and selection

Primary structures must be experimentally determined protein-DNA complexes
(X-ray, NMR, or cryo-EM) with a directly modelled TF-DNA interface. AlphaFold,
AF3, homology models, and DNA-free structures are not eligible for the primary
challenge. Candidate PDBs are selected from official RCSB/wwPDB records using
TF identity and DNA-bound annotations, before reading model-vs-label results.

For a TF with multiple candidates, the deterministic priority is:

1. correct TF/construct and a direct DNA-bound complex;
2. complete DNA-binding domain and protein-DNA interface;
3. an unambiguous sequence-compatible DNA duplex;
4. experimental resolution/method quality;
5. higher interface completeness.

The primary PDB is fixed by the structure audit. Alternative PDBs can be
listed as pre-specified sensitivity analyses only. No PDB, chain, orientation,
or register may be selected using any SaMBA effect, correlation, rank, or
pairwise result.

## 3. DNA mapping

Mapping uses only PDB polymer sequences, residue/chain metadata, exact
substring alignment, reverse-complement consistency, and the deterministic
canonical motif/register rule below. A SaMBA site is eligible only when its
canonical mutation positions map to a unique PDB duplex position set. Multiple
valid registers, missing duplex strands, or unresolved orientation are
`MAPPING_AMBIGUOUS` and excluded from the primary challenge. Register search is
never optimized against labels or predictions.

The canonical rule is: identify the SaMBA reference sequence and its mutated
base-pair coordinates, derive the contiguous canonical site motif from those
coordinates, and compare it to each PDB DNA duplex strand and its reverse
complement. An exact substring is preferred. If no exact substring exists, a
label-blind global sequence alignment is used with match `+2`, mismatch `-1`,
and gap `-2`; the orientation with the higher alignment score is retained and
ties are `MAPPING_AMBIGUOUS`. A mapping is retained only when every mutated
base-pair coordinate maps to a nucleotide in the PDB strand, the aligned
Watson-Crick partner is present, and the best alignment is unique. Flanking
SaMBA sequence outside the PDB duplex is not scored. The alignment algorithm
and tie rule are fixed before model inference and never use effects or
predictions.

## 4. DeepPBS method and score

Primary predictions use the repository's existing official DeepPBS checkout at
commit `8bfb211dd67f02877841f6f33aa493ddf7daedf9`, the existing five-checkpoint
ensemble (`828`, `529`, `173`, `898`, `820`), the checked-in preprocessing and
prediction configuration, and the Linux CPU environment documented in
`external/deeppbs/RUNBOOK.md` and `V1_3_HOLDOUT_PREDICTION_PROTOCOL.md`.

No checkpoint, model version, preprocessing, sign, or implementation may be
changed after seeing natural-TF performance. If the frozen environment cannot
reproduce an eligible structure, the site is `NOT_EVALUABLE`; a replacement
DeepPBS version is secondary/post-hoc only.

For each mapped position and each of the three non-wild-type base pairs, the
prediction is the frozen native PWM delta-logP score above. The model output is
not treated as Kd or calibrated affinity.

## 5. Primary endpoints

Metrics are exactly those in `V1_3_EVALUATION_CONTRACT.md`:

- per TF/site global mutation-effect Spearman (Pearson descriptive only when
  finite and meaningful);
- position sensitivity using `mean_abs_effect` and its Spearman correlation;
- within-position residual Spearman after subtracting the position mean from
  experimental and predicted effects;
- within-position pairwise accuracy, with the frozen tie handling and chance
  baseline 0.5.

The statistical unit is TF/site. Report per-site values, then first aggregate
sites within TF, and finally summarize across TFs using median and bootstrap
intervals where feasible. Mutation rows are not independent biological
replicates.

## 6. Two-stage boundary

Stage A writes only `results/v1_3/natural_tf_predictions_unscored.tsv` and
`reports/v1.3/V1_3_3_PREDICTION_FREEZE.md`. It may read structure metadata and
SaMBA sequence/mapping metadata, but it must not load SaMBA effect values or
compute any prediction-vs-label metric. The structure/mapping audit and Stage A
prediction are committed before reveal.

Stage B may join frozen predictions to the canonical labels and write
`natural_tf_decomposition*.tsv`, evaluation reports, and figures. No endpoint,
eligibility rule, PDB, register, orientation, or score sign may be changed
after reveal.

## 7. Interpretation rules frozen in advance

- Natural identity high while designed identity is near chance: support a
  designed-protein domain-shift interpretation.
- SaMBA technical identity repeatability high but natural DeepPBS identity low,
  with designed identity low: support a structure-model identity limitation;
  this is the gate for considering v1.4 identity-aware modelling.
- Natural global/position/identity all poor: do not call this identity-specific
  failure; retain assay, structure/site mapping, and domain-transfer alternatives.
- Too few eligible structures or unambiguous mappings: report
  `STRUCTURAL_COVERAGE_LIMITED` and keep v1.4 closed.

SAMPDI-3D/3Dv2 and motif baselines are secondary. SaMBA technical spots are
technical repeatability, not biological replicates. No natural-TF result can
retroactively alter the frozen designed benchmark.
