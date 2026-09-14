# Structure-to-7-mer Scoring Protocol (Frozen)

This protocol is frozen before reading any v1.0 PBM performance. It is applied identically to NA-MPNN PPMs and, where available, DeepPBS PWM outputs.

1. The candidate universe is the 8,192 canonical reverse-complement classes defined by the frozen v0.6 sequence-equivalence utility.
2. A PPM/PWM is converted to a log-probability score using `sum(log(P(base_i)+1e-12))`.
3. For a structure with multiple contiguous DNA runs, all structurally valid contiguous 7-bp windows and both candidate orientations are scored; the maximum is retained. This is a fixed register aggregation rule, not a per-protein PBM-driven search.
4. The experimental/design core is taken from the frozen target metadata (`GCAGG`, `GCAGGA`, `TGCACA`, or `CTGACG`). No register is selected by maximizing Spearman.
5. Missing residues, ambiguous chain mapping, or a DNA run shorter than 7 bp produce `not evaluable`, never an imputed score.
6. PPM values are ranking proxies; they are not numerically equated with PBM E-scores.

The output records the selected structural window and orientation so every ranking is auditable.
