# DeepPBS v0.4.2 Scoring Protocol

This protocol is frozen before any designed-DBP DeepPBS result is evaluated.

The upstream semantics supporting this protocol are documented in
`docs/v0_4_2/DEEPPBS_OUTPUT_SEMANTICS.md`, with source-level references to
`deeppbs/dna_encodings.py`, `run/predict.py`, and
`run/models/model_v2.py` at DeepPBS commit
`8bfb211dd67f02877841f6f33aa493ddf7daedf9`.

## Inputs

Use only public/original protein-DNA complex structures listed in `metadata/v0_4_2/designed_structure_manifest_v2.csv`.

## Primary mapping

If official DeepPBS returns a position-wise PWM/profile over a structure-defined DNA binding window, the primary candidate 7-mer score is:

`score(D) = max_offset sum_i log(P(D_i at aligned position offset+i) + 1e-9)`

The offset range is all contiguous 7-bp windows fully contained in the structure-defined DeepPBS output window. Reverse-complement candidate classes are scored by the maximum over both orientations.

## Constraints

- Do not tune offsets, checkpoint, or aggregation after inspecting PBM E-scores.
- Do not call the PWM-derived value a PBM affinity, Kd, binding probability, or
  calibrated binding score.
- If fewer than 4 designed DBPs are evaluable, DeepPBS remains diagnostic.
