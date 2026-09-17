# v1.3.2 External Validation Contract

Status: **FROZEN 2026-09-17, before external data download or result inspection**.

This contract extends, but does not modify, `V1_3_EVALUATION_CONTRACT.md`. The designed-DBP development and locked-holdout results remain frozen. External results cannot change the primary metric definitions, perturbation grouping, inclusion rules, or method sign. Changes are permitted only for a demonstrated implementation or source-interpretation bug and must be recorded as an amendment before rerunning an affected endpoint.

## Questions and claim boundary

External Dataset A (Afek et al. SaMBA, DOI `10.1038/s41586-020-2843-2`) asks whether identity-level experimental effects have sufficient **technical spot-level repeatability**. Array spots are technical observations, not biological replicates or an independent experimental noise ceiling.

External Dataset B (SAMPDI-3D T227, DOI `10.1093/bioinformatics/btab567`) asks whether global, position, and within-position identity model performance can be evaluated in natural protein-DNA complexes. It is not assumed in advance that T227 contains multiple mutant identities per position.

No external dataset will be described as validating a universal `position >> identity` phenomenon unless its frozen endpoints are evaluable and its results support that statement. Negative or non-evaluable results are retained.

## Source and provenance rules

Publisher supplements, NCBI GEO records/files, official method distributions, and author repositories are preferred in that order. Figure OCR and manual transcription from plotted points are prohibited. New source files are stored below `data/raw/v1_3_2_external/`.

Each newly downloaded source file is hashed once through the existing `scripts/v1_3/file_manifest.py` cache. Later runs compare path, size, and mtime and reuse the cached SHA256 when unchanged. Generated tables are not repeatedly fingerprinted, and no repository-wide hash is created.

## Dataset A: perturbation strata

Every SaMBA record must be classified from source metadata as one of:

- `canonical_watson_crick`: replacement of one canonical Watson-Crick base pair by another canonical Watson-Crick base pair.
- `mismatch`: a non-Watson-Crick paired perturbation.
- `other_or_unresolved`: chemistry or pairing cannot be assigned without inference.

The primary SaMBA analysis uses only `canonical_watson_crick`. Mismatches are a separate secondary analysis. The two strata are never pooled to increase sample size. `other_or_unresolved` is audited but excluded from both analyses.

The mutation key is `(tf_id, binding_site_id, position, wt_basepair, perturbation)`. Strand/orientation is taken from documented source metadata. Reverse-complement or base-pair canonicalization may standardize notation but may not collapse biologically distinct mismatches.

## Dataset A: spot-level table and signal

The processed spot table retains one row per recoverable spot with at least:

`tf_id`, `binding_site_id`, `position`, `wt_basepair`, `perturbation`, `perturbation_type`, `spot_id`, `array_id`, `raw_signal`, `normalized_signal`, and `source_file`.

Official spot-level normalized signal is primary when its definition is documented. If only spot intensity is supplied, each split-specific perturbation estimate is the log2 ratio of its mean valid intensity to the mean valid intensity of the matched canonical reference sequence on the same array or documented normalization block. No value is imputed. Non-positive intensity, absent reference, saturation flag, or unresolved normalization produces an explicit exclusion reason. If the source supplies only sequence-level processed values, spot-split repeatability is `NOT_EVALUABLE` rather than reverse-engineered.

## Dataset A: frozen technical split

The primary split seed is `1301`. Rows and grouping keys are sorted lexicographically before assignment. Within each `(tf_id, binding_site_id, sequence_or_perturbation_key, array_id)` stratum, a NumPy `PCG64` permutation is drawn from one generator initialized once with seed 1301; permuted observations are assigned alternately to Split A and Split B. A stratum requires at least two valid spots. When `array_id` is genuinely unavailable, all spots for the sequence form one explicitly labeled unknown-array stratum. Each split must retain at least one spot for every included perturbation and its matched reference.

Seeds 1301 through 1400 are reserved for a secondary repeated-split stability analysis. They cannot replace the primary split. Spot IDs, array IDs, or source row IDs used for assignment are preserved so the split is reproducible.

## Frozen decomposition metrics

For both external datasets, larger effect values are oriented toward stronger binding only when the source definition supports that orientation. If the published endpoint is a destabilizing energy, both experimental and prediction values may be retained in that documented direction; any sign transform must be source-defined and recorded before method correlations are inspected.

For an effect matrix `y(p,b)`:

1. **Global:** Spearman and Pearson across mutation rows within TF/site or complex.
2. **Position:** `mean_abs_effect(p) = mean_b(abs(y(p,b)))`, followed by Spearman and Pearson across positions.
3. **Identity:** `delta(p,b) = y(p,b) - mean_b(y(p,b))`, followed by pooled residual Spearman and Pearson within TF/site or complex.
4. **Identity pairwise:** within-position pairwise accuracy, with 1 for a correct strict order, 0 for an incorrect strict order, and 0.5 for a prediction/estimate tie. Experimental ties are excluded. Chance is 0.5. Kendall-style concordance is `2 * accuracy - 1`.

For SaMBA, `y_A` and `y_B` are the two independently aggregated technical spot estimates, so the metrics quantify technical repeatability. For T227, the same metrics compare a published prediction with experiment when both are available.

Position-label and within-position identity-label permutation use the v1.3 definitions, 2,000 permutations, and seed 1301. Mutation rows are not treated as independent TF replicates.

## Frozen evaluability thresholds

Within each TF/site or complex:

- Global correlation requires at least 6 complete perturbation rows.
- Position correlation requires at least 3 eligible positions.
- Identity residual correlation requires at least 6 complete rows distributed over at least 3 positions, with at least 2 mutant identities at every included position.
- Identity pairwise accuracy uses the same eligible identity positions and requires at least 3 non-tied experimental pairs in total.

Rows outside a metric's eligible set are not borrowed from mismatch/other strata. A dataset may be globally evaluable but identity-level `NOT_EVALUABLE`.

SaMBA site-level metrics are first summarized by the median within TF when a TF has multiple sites. Across-TF summaries then report per-TF values, `n_tf`, median, and mean. T227 metrics are calculated per complex; multiple complexes from the same TF are summarized by their within-TF median before across-TF summaries. Technical spots and mutation rows never inflate the biological unit count.

## Landscape variance decomposition

For each sufficiently complete effect matrix, define the overall mean `mu`, position component `alpha_p = mean_b(y(p,b)) - mu`, and residual `epsilon(p,b) = y(p,b) - mu - alpha_p`. Using population variance (`ddof=0`) and row-expanding `alpha_p` to the observed mutation rows:

`fraction_position = Var(alpha_p expanded) / Var(y)`

`fraction_identity = Var(epsilon) / Var(y)`

These describe experimental or predicted **landscape variance allocation**, not predictive performance and not the percentage of correlation explained by position. They are reported only when total variance is positive and the identity evaluability matrix is satisfied. No subtraction of Spearman correlations is interpreted as explained signal.

## T227 method handling

Published SAMPDI-3D and FoldX prediction values are preferred over rerunning software. A method is included only when prediction, experimental effect, complex, position, wild-type base pair, and mutant base pair can be joined without outcome-driven mapping. No sign, register, chain, subset, or mutation grouping is selected by correlation.

If published predictions cannot be recovered, the dataset audit and experimental structure analysis proceed and the method is marked `NOT_REPRODUCED` in `reports/v1.3/METHOD_REPRODUCIBILITY.md`. Environment troubleshooting is bounded and cannot delay dataset-level conclusions.

## Uncertainty and reporting

Per-TF/site or per-complex values are always retained. Across-dataset summaries use TF as the primary unit, with deterministic TF-level bootstrap intervals only when at least five TFs are evaluable; otherwise intervals are omitted and the small `n` is explicit. Technical split stability is secondary and cannot be described as biological reproducibility.

Primary canonical SaMBA, secondary mismatch SaMBA, and T227 outputs remain separate. No TF, site, position, or method is excluded after viewing performance except for a pre-specified missingness/evaluability rule recorded with a reason.

## v1.4 gate

The v1.4 model gate opens only if external technical identity repeatability is demonstrable, natural mutation data contain evaluable identity structure, and frozen model identity prediction remains weak. If technical identity repeatability is weak, the gate stays closed and the interpretation shifts toward benchmark identifiability. If the phenomenon is confined to designed DBPs, claims are narrowed accordingly. No model is trained before this gate decision.
