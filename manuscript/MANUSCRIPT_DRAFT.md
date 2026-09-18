# Disentangling Protein Conditionality, Position Sensitivity, and Nucleotide Identity in DNA-Binding Specificity Prediction

## Abstract

DNA-binding specificity predictors are commonly summarized with one aggregate correlation across all mutations. That summary is useful, but it can combine two distinct questions: whether a model identifies sensitive positions and whether it identifies which nucleotide is preferred at each position. We evaluated these questions in a frozen designed-DNA-binding-protein (DBP) benchmark using global ranking, position-level ranking, within-position residual ranking, and identity pairwise accuracy. We also tested protein conditionality with protein- and target-shuffle diagnostics, evaluated a prospective locked holdout, decomposed experimental and predicted landscape variance, and constructed a position-only counterfactual that removes within-position nucleotide identity from predictions. DeepPBS showed substantially stronger position-level than identity-level agreement in both development proteins and the three locked holdouts. Across development proteins, the median global, position, and identity-residual Spearman values were approximately 0.354, 0.406, and 0.114; the locked holdout medians were 0.085, 0.560, and 0.001, with pairwise accuracy 0.548. Removing within-position identity from predictions sometimes preserved or increased aggregate correlation, demonstrating that aggregate mutation-effect performance can remain similar after identity discrimination is removed. Experimental landscapes contained nonzero within-position variance, while DeepPBS often allocated a larger share of prediction variance to between-position effects. Independent SaMBA technical spot splits showed strong identity-level repeatability (median residual Spearman 0.823 and pairwise accuracy 0.814), establishing that nucleotide identity can be technically repeatable in a distinct assay without defining a competition-assay noise ceiling. These results motivate benchmarks that report protein conditionality, position sensitivity, and nucleotide identity as separate axes.

## Introduction

Predicting protein-DNA specificity is central to interpreting regulatory variants, designing DNA-binding proteins, and estimating off-target binding. Modern predictors are often evaluated by correlating predicted and measured effects across a mutation landscape. A single global correlation is attractive because it is compact and comparable across proteins, but it does not identify the source of agreement. A model can rank positions that are generally sensitive while assigning similar effects to all four mutant bases at each position. Such a model may obtain a respectable aggregate score while failing the mechanistic question that distinguishes alternative nucleotides.

The distinction matters for designed DBPs. The benchmark used here contains matched protein and competition-landscape measurements, allowing the evaluation to separate protein conditionality from DNA mutation resolution. Protein conditionality asks whether predictions change when the protein or target context changes. DNA mutation resolution asks whether the model captures global effects, position sensitivity, and nucleotide-specific residual effects. These are orthogonal axes: a model can be protein-conditioned but position-dominant, or position-sensitive but weakly identity-aware.

We therefore froze a metric contract before reveal of a three-protein prospective holdout and applied the same contract to the development set. Global Spearman correlation measures aggregate ranking. Position Spearman correlates position-averaged effects. Identity residual Spearman compares deviations from each position mean, and identity pairwise accuracy asks whether the model orders two mutant bases at the same position correctly. We complement these predictive metrics with a variance allocation that describes the experimental or predicted landscape itself, not the fraction of predictive performance attributable to a feature. Finally, a position-only counterfactual replaces each predicted mutation effect with the mean prediction at its position, preserving position information while deleting nucleotide identity.

## Results

### Protein-conditioned models exhibit conditionality collapse

Cached M1/M2/M3 diagnostics were reused without retraining. Protein-shuffle and target-shuffle perturbations test whether predictions depend on the intended protein and target context. The M3 conditionality diagnostics show near-invariant predictions under protein replacement in the cached developmental analysis, with protein-shuffle correlations close to one and a negligible conditioning effect. This is a conditionality failure mode distinct from the position-versus-identity failure described below. The result is reported as a diagnostic of the existing model family, not as evidence that all protein-DNA predictors lack conditionality.

### DeepPBS is more position-sensitive than identity-sensitive

The designed development set gives a median global Spearman of approximately 0.354 and a median position Spearman of approximately 0.406. The median identity-residual Spearman is approximately 0.114 and identity pairwise accuracy is approximately 0.506. Thus position-level agreement is stronger than within-position nucleotide discrimination. The per-protein table is in `results/v1_3/decomposition_per_protein.tsv`; the publication-oriented source table is `results/v1_3/figure2_revised_data.tsv`.

The variance analysis describes the corresponding landscapes. For every protein, we calculated the total variance of the mutation matrix, the variance of the position component, and the variance of the residual after subtracting the position mean. The resulting fractions are reported separately for experimental competition landscapes and DeepPBS predictions in `results/v1_3/variance_decomposition_per_protein.tsv`. This is a descriptive allocation of landscape variance. It is not a claim that a given percentage of model performance is caused by position.

### Prospective locked holdout reproduces position >> identity

DBP023, DBP056, and DBP062 were predicted before their labels were joined. Their median global Spearman is approximately 0.085, median position Spearman is approximately 0.560, and median identity-residual Spearman is approximately 0.001. Median identity pairwise accuracy is approximately 0.548, only modestly above the fixed 0.5 chance baseline. All three proteins therefore directionally reproduce the development result: position sensitivity is stronger than nucleotide-identity discrimination. The holdout is small and supports a directional replication, not a broad significance claim.

### Position-only counterfactuals can preserve aggregate performance

For each protein, we formed `y_position_only(p,b) = mean_b y_hat(p,b)`. This counterfactual preserves the predicted position effect and removes all within-position mutant-base differences. It was compared with the experimental landscape using the same global Spearman metric. Results are in `results/v1_3/position_only_counterfactual.tsv` and Figure 3.

The counterfactual does not uniformly improve or worsen the score. For DBP001, DBP006, DBP009, and DBP023 it is at least as high as the full score; for other proteins it decreases. The key result is therefore not that the global metric is useless, nor that all performance comes from position. Rather, aggregate mutation-effect performance can remain similar after removing within-position nucleotide-identity discrimination. A global score can consequently mask a mechanistically important identity-level failure.

### DeepPBS often allocates more predicted variance to positions than identity

Experimental landscapes contain substantial between-position variance, but they also contain a measurable within-position residual component. DeepPBS predicted landscapes commonly allocate an even larger fraction to the position component, especially in the prospective holdout. The per-protein values, including total, position, and residual variances, are provided in the variance-decomposition table. This comparison addresses the question “does the model allocate too much predicted variance to between-position effects and too little to identity?” without converting variance fractions into claims about predictive performance.

### SaMBA measurements show independent technical identity repeatability

The canonical SaMBA technical spot split provides an independent assay context. Across seven TF-level canonical analyses, median global Spearman is 0.947, position Spearman 0.881, identity-residual Spearman 0.823, and identity pairwise accuracy 0.814. These are technical repeatability results from spot splits, not biological replicates and not a competition-assay noise ceiling. SaMBA and the competition benchmark estimate different observables; the result only establishes that identity-level effects can be technically repeatable in another assay.

### Natural-TF structural challenge remains not evaluable

The frozen PDB 2STT input is a 25-model NMR ensemble with 2,670 atoms per model, protein chain C, DNA chains A and B, no alternate locations, and no water or heteroatom records. The historical failure occurred during DeepPBS preprocessing before an NPZ was produced. The deterministic MODEL 1 extraction is stored at `data/processed/v1_3_4/2STT_model1.pdb` and is recorded in the forensic audit. The positive control passed using an existing DBP035 NPZ. The bounded local rescue then failed at import because the current environment lacks `torch_cluster`; historical Docker/WSL logs already document unavailable container execution. No alternate structure, conformer cherry-picking, webserver result, or new model was used. The natural-TF primary result therefore remains `NOT_EVALUABLE`. Any future successful MODEL 1 run would be labelled `POST-REVEAL TECHNICAL SENSITIVITY`, not prospective validation.

## Discussion

The central implication is methodological. Aggregate correlations can reward a model for identifying where a landscape changes while concealing whether it knows which nucleotide causes the change. A position-only counterfactual makes this failure mode explicit without inventing a new primary metric. Reporting global, position, residual identity, and pairwise results together makes it possible to distinguish broad ranking from nucleotide-level resolution.

The prospective holdout strengthens the position-dominant observation because its predictions were frozen before reveal. It does not eliminate uncertainty: the holdout contains only three proteins, and the competition source provides published means rather than replicate-resolved measurements. We therefore cannot assign the weak identity result uniquely to model error versus competition-assay identifiability. SaMBA technical repeatability is informative but not a substitute for a competition noise ceiling.

Protein conditionality is a second, orthogonal axis. A predictor can resolve positions while failing to condition on protein identity, and it can condition on protein identity while remaining identity-blind at a position. The M1/M2/M3 shuffle diagnostics should therefore remain alongside the mutation-resolution decomposition in future benchmarks.

The DBP35/DBP35opt comparison remains a context-confounded case study because competitor and target concentrations differ. Landscape similarity can be described, but protein mutations cannot be assigned a causal effect without matched assay context. The natural-TF challenge also remains unresolved after the bounded rescue. Consequently, the current evidence supports a carefully scoped designed-DBP conclusion rather than a universal natural-to-designed generalization claim.

## Methods

### Designed DBP data and split contract

The analysis used the official competition means parsed into `data/processed/v1_3_competition_mutations_official_mean.parquet`. Seven DBPs were treated as development proteins and DBP023, DBP056, and DBP062 as a prospective locked holdout. Holdout DeepPBS predictions were generated and frozen before labels were joined. No new neural network, GNN, fine-tuning, or hyperparameter search was performed in v1.3.4.

### Prediction metrics

Global Spearman is the rank correlation across all mutation rows for a protein. Position Spearman correlates position-averaged experimental and predicted effects. For identity analysis, each effect is residualized by subtracting the mean across mutant bases at the same position. Identity residual Spearman is the rank correlation of these residuals. Identity pairwise accuracy is the fraction of within-position mutant-base pairs ordered concordantly, with ties scored 0.5. The fixed pairwise chance baseline is 0.5. Permutation and bootstrap units follow the frozen v1.3 evaluation contract and are not re-optimized here.

### Landscape variance decomposition

For a landscape `y(p,b)`, `mu` is the grand mean, `alpha_p = mean_b(y(p,b)) - mu`, and `epsilon(p,b) = y(p,b) - mu - alpha_p`. We report total variance, variance of `alpha`, variance of `epsilon`, and each component divided by total variance. This analysis is performed separately for experimental competition landscapes and DeepPBS predictions, per protein and cohort.

### Position-only counterfactual

For each predicted mutation, the position-only value is the mean predicted effect across mutant bases at that position. Experimental values are left unchanged. Global Spearman is recomputed against the experimental landscape. The full identity residual Spearman and pairwise accuracy are retained as diagnostics; no new primary endpoint is introduced.

### SaMBA processing

SaMBA values are the canonical Watson-Crick technical spot split with the frozen seed and source-integrity exclusions documented in the v1.3.2 audit. Spot splits are technical repeatability analyses, not biological replicate analyses.

### Structural audit and provenance

The 2STT forensic audit records file size, model and atom counts, chains, residue names, alternate locations, waters, heteroatoms, formatting, logs, process exit code, expected NPZ path, and actual output. MODEL 1 was extracted deterministically without minimization, remodeling, mutation, fitting, or label-guided register changes. File hashes use the persistent size/mtime-aware manifest; unchanged files are not rehashed.

## Data and code availability

All result tables, figure source files, audit reports, and scripts are retained in this repository. The main analysis outputs are `results/v1_3/variance_decomposition_per_protein.tsv`, `results/v1_3/position_only_counterfactual.tsv`, and `results/v1_3/natural_tf_2stt_rescue_status.tsv`.

## Conclusion

DeepPBS reproduces position-level sensitivity more reliably than within-position nucleotide identity in the frozen designed-DBP benchmark and its prospective holdout. Position-only counterfactuals show why aggregate scores alone are insufficient. Independent SaMBA data demonstrate that identity can be technically repeatable, but do not provide a competition noise ceiling. The evidence is sufficient for a consolidated, scope-limited paper and for planning an identity-aware model, but the v1.4 training gate remains closed until a broader defensible training set and replicate-resolved competition evidence are available.

## Figure captions

**Figure 1.** Conceptual axes separating protein conditionality from DNA mutation resolution.

**Figure 2.** Designed DBP development and prospective locked holdout metrics. Holdout predictions were frozen before reveal.

**Figure 3.** Aggregate metrics can mask loss of nucleotide identity. Full versus position-only global scores, experimental and predicted position-variance fractions, and identity diagnostics are shown.

**Figure 4.** SaMBA technical repeatability for global, position, identity-residual, and pairwise endpoints.

**Figure 5.** Cached protein-conditionality diagnostics for M1, M2, and M3; no retraining was performed.
