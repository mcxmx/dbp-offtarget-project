# v1.3.4 DeepPBS Failure Audit

Audit status: completed from the frozen input and existing logs; no labels or
experimental performance were used to choose a structure or debug route.

## Frozen input

| field | observed value |
|---|---|
| input | data/raw/v1_3_3_external/structures/2STT.pdb |
| file size | 5,446,764 bytes |
| MODEL records | 25 |
| ENDMDL records | 25 |
| conformers | 25 (MODEL 1 through MODEL 25) |
| atoms per MODEL | 2,670 for every model |
| MODEL 1 chains | A, B, C |
| DNA chains | A, B |
| protein chain | C (frozen mapping) |
| nucleotide residue names | DA, DC, DG, DT |
| alternate locations | none observed; blank altloc only |
| waters / HETATM in MODEL 1 | none |
| total ATOM/HETATM records | 66,750 |
| unusual formatting | fixed-width ATOM records; no malformed records observed in the audited atom block |

The file is therefore confirmed as a multi-model NMR ensemble, not a single
structure. The deterministic rescue uses MODEL 1 and does not cherry-pick a
conformer by SaMBA performance.

## Existing failure evidence

The v1.3.3 freeze records the official DeepPBS checkout at commit
8bfb211dd67f02877841f6f33aa493ddf7daedf9, bundled checkpoints
828,529,173,898,820, and the fixed full 2STT input. The documented local
preprocessing returned without producing an NPZ; the frozen prediction rows
therefore carry NOT_REPRODUCED_DEEPPBS_PREPROCESSING_NO_NPZ. The recorded
artifact path was results/v1_3/empty_npz/2STT_Ets1_site_1.npz_predict.npz
and analogous site paths; no valid NPZ was generated.

The older environment log records Docker unavailable (exit 127), WSL without
an installed distribution (exit 50), and a Windows process launch error. Those
facts identify an environment/implementation reproducibility boundary, but do
not by themselves show whether the 2STT input is intrinsically unsupported.

## Failure-stage conclusion

The primary failure stage is DeepPBS preprocessing / NPZ generation, before
prediction and before any metric calculation. The failure is technical and the
natural-TF primary challenge remains NOT_EVALUABLE; it is not a negative model
result. A successful MODEL 1 rescue would be post-reveal technical sensitivity
only.

## Provenance boundaries

The frozen 2STT structure is not replaced by 1HLO, 1R4R, 1TGH, another ETS1
structure, AlphaFold, or any predicted structure. No all-conformer search,
register search, orientation change, checkpoint change, sign change, or
SaMBA-guided debugging is permitted.
