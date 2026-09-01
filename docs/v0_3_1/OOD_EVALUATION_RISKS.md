# OOD Evaluation Risks

Audit date: 2026-09-01

## Core Risk

If future models are trained on natural DBP data from one assay type and evaluated on GSE237017 designed DBPs from uPBM, performance drop may reflect both:

- protein distribution shift
- assay distribution shift

Therefore, natural-to-designed performance drop alone is not sufficient evidence that designed DBPs are intrinsically harder or uniquely OOD.

## Concrete Confounding Example

Natural training data:

- HT-SELEX enrichment

Designed external test:

- uPBM processed E-score

Observed degradation could be caused by differences in sequence library, readout, score scale, preprocessing, dynamic range, or noise, not only designed-protein biology.

## Required Controls for v0.4/v0.5

Preferred benchmark layout:

- Natural PBM/uPBM train to natural PBM/uPBM held-out.
- Natural PBM/uPBM train to designed uPBM external.
- Natural HT-SELEX train to natural HT-SELEX held-out.
- Optional cross-assay analysis reported separately.

## Interpretation Rule

Do not claim a designed-DBP OOD failure unless an assay-matched natural control shows that the drop exceeds ordinary held-out natural PBM/uPBM performance loss.

## Current Status

v0.3.1 does not download a large natural PBM dataset. It only records this risk and keeps GSE237017 positioned as a conditional external benchmark.
