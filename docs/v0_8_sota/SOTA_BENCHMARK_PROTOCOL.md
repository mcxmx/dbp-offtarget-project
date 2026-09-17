# V0.8 SOTA Benchmark Protocol

This is a frozen audit of already-run structure diagnostics plus historical sequence baselines. No model is fine-tuned on GSE237017, no hyperparameter is selected from the exposed cohort, and no new sequence architecture is introduced.

## Primary unit

The independent statistical unit is protein (N=7 for the designed cohort). Candidate-level scores (8,192 canonical RC classes per protein) are used only to calculate within-protein Spearman and ranking metrics. Structure methods are reported with their actual coverage (N=2 for DeepPBS and NA-MPNN), never as a seven-protein estimate.

## Uniform scoring

DeepPBS PWM and NA-MPNN PPM outputs are converted to a fixed canonical 7-mer log-probability ranking. No PBM E-score numerical equivalence is claimed. Register/window and strand conventions are inherited from the frozen v0.4 artifacts. Missing structures are missing, not imputed.

## Conditionality

Protein/target shuffle diagnostics are reported only where the model defines a counterfactual. Structure scores are not assigned meaningless shuffle values because each structure is protein-specific. A high Spearman with no conditionality evidence is classified as collapsed, not successful specificity learning.

## Exposure

All seven GSE237017 DBPs are development-exposed. These results are benchmark/development diagnostics and cannot be called external validation.
