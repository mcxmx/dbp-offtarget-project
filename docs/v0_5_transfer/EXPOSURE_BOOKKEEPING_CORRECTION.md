# v0.5 Designed-DBP Exposure Correction

## Current status

The complete v0.5 primary evaluation inspected all seven designed DBPs.
Therefore no GSE237017 designed protein is currently an untouched confirmatory
test protein.

The earlier phase records remain valid historical snapshots:

- v0.5 fold-1 engineering smoke exposed DBP1 and DBP3.
- Phase 5A local smoke exposed DBP5 and DBP35.
- Phase 6A dense smoke exposed DBP48.
- The complete v0.5 primary evaluation subsequently exposed DBP6 and DBP9 as
  well as re-exposing the other five proteins.

The historical `phase6a_exposure_status` fields are intentionally not
rewritten. Current status is recorded in
`metadata/v0_5_transfer/exposure_manifest.csv`.

## Consequence

All later GSE237017 analyses are developmental/exploratory. A future final
model claim requires an independent designed-DBP dataset or a prospective
holdout that was not inspected during model development.
