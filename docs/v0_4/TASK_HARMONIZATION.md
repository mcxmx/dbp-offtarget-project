# v0.4 Task Harmonization

The ground truth task is per-protein ranking of reverse-complement canonical 7-mers by processed experimental uPBM E-score consensus.

## PBM E-score

`experimental_escore_consensus` is a processed experimental uPBM E-score consensus. It is suitable here as a per-protein ranking target. It is not Kd, binding free energy, binding probability, or an absolute cross-protein affinity scale.

## Sequence-Only Baselines

Hamming/edit/k-mer metrics are sequence-only proxy metrics computed against the v0.3.1 paper motif reference. They are not protein-conditioned.

## DeepPBS

DeepPBS predicts a structure-conditioned DNA PWM. A fair mapping to this benchmark would require a validated protein-DNA complex structure for each designed DBP and a predeclared PWM-to-7mer scoring rule. v0.4 does not generate DeepPBS predictions because the official preprocessing chain is not runnable in the current environment.

## NA-MPNN

NA-MPNN specificity mode predicts a PPM over residue types at nucleic-acid positions in a supplied structure. For DBP35 and DBP48, v0.4 maps the predicted DNA-position probabilities to canonical 7-mer scores as follows:

1. Use the official `s_70114.pt` specificity checkpoint.
2. Extract predicted probabilities for DA/DC/DG/DT at nucleic-acid positions from the official inference `.npz` output.
3. Preserve each contiguous DNA-chain run rather than concatenating unrelated chains.
4. For each candidate 7-mer, compute the sum of log probabilities in every complete 7-position window.
5. Score a canonical RC class by the maximum score over candidate orientation and reverse-complement orientation.

The resulting score is named `partial_structural_ppm_best_window_log_probability`. It is a derived ranking score, not a PBM E-score and not an affinity. DBP48/8TAC has NA-MPNN validation-split overlap and must be interpreted as a diagnostic, not a zero-shot result.

## Main Metrics

Metrics are computed per protein and macro-summarized. The primary ranking metric is Spearman correlation; NDCG@1%, NDCG@5%, top-1% recovery, and sampled pairwise ranking accuracy are secondary. Rows from different proteins are not pooled into a single main correlation.
