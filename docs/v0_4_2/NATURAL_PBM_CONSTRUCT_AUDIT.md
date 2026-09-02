# v0.4.2 Natural PBM Construct Audit

Audit date: 2026-09-02

The v0.4.1 natural PBM benchmark recovered protein sequences primarily from UniProt/reference records. The local metadata explicitly does not claim those sequences are the exact PBM assay constructs.

## Result

- Natural proteins audited: 57
- Exact experimental constructs recovered: 0
- Domain/truncated constructs reconstructed from reported coordinates: 0
- Confirmed full-length assay constructs: 0
- Unknown assay constructs: 57
- High-confidence construct coverage: 0.000

## Consequence

`FULL_LENGTH_REFERENCE` remains a sensitivity benchmark because it is reproducible and has sequence provenance, but it is not an assay-aligned construct benchmark. `ASSAY_ALIGNED_PROTEIN` is empty in v0.4.2 because no construct sequence or reliable construct coordinate provenance was recovered.

No missing construct sequence is filled by guessing, domain heuristics, or family-level substitution.
