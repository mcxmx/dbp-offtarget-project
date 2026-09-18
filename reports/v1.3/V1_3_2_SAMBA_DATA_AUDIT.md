# v1.3.2 SaMBA Data Audit

Status: **complete; source and mapping audit frozen 2026-09-17**.

## Source inventory

The official Afek et al. supplementary files were downloaded from Springer Nature and NCBI GEO:

| source | local file | role |
|---|---|---|
| Supplementary Table 1, DOI `10.1038/s41586-020-2843-2` | `data/raw/v1_3_2_external/samba/41586_2020_2843_MOESM4_ESM.xlsx` | binding-site inventory, SaMBA processed sheets, and raw mutation-array signal columns |
| Supplementary Table 4 | `data/raw/v1_3_2_external/samba/41586_2020_2843_MOESM7_ESM.xlsx` | 12-site calibration workbook with canonical/mismatch site metadata and raw spot-level signal |
| GEO `GSE156375` raw archive | `data/raw/v1_3_2_external/samba/GSE156375_RAW.tar` | 29 processed array sample files; retained for provenance and follow-up parser work |

The persistent SHA256 entries and publisher/NCBI URLs are recorded in `reports/v1.3/V1_3_2_EXTERNAL_DATA_PROVENANCE.md` and `artifacts/cache/file_manifest.json`.

The paper and GEO metadata describe approximately 8-20 technical spots per sequence. These are treated as technical within-array/assay observations, not biological replicates.

## Parsed calibration inventory

The deterministic parser `scripts/v1_3/run_v1_3_2_external.py` recovered 22,026 raw signal rows from Supplementary Table 4. Before source-integrity filtering, the workbook contains 12 sites from 7 TFs (`Ets1`, `TBP`, `Max`, `Cbf1`, `Egr1`, `p53`, `GR`), 2,125 distinct sequences, 414 canonical Watson-Crick mutant sequences, 1,683 mismatch sequences, and 12 reference sequences. The processed spot table is `data/processed/v1_3_2_samba_spots.parquet`.

Canonical and mismatch perturbations are retained in separate strata. A canonical mutation changes one Watson-Crick base pair to another Watson-Crick base pair. A mismatch is a non-Watson-Crick paired perturbation. No mismatch row is pooled into the primary canonical analysis.

## Source-integrity and mapping audit

The same 60-mer sequence occurs in multiple named probe blocks for several sites. The parser therefore does not choose a block by maximizing agreement with a result. It instead reconstructs the workbook's published median from every raw row joined by sequence. A site is included only when every sequence's raw median exactly equals the published `Median_signal` (within floating-point tolerance).

Ten sites pass this check. Ets1 site 2 and Ets1 site 4 contain sequence-level raw/published median inconsistencies and are retained in the spot table with `source_mapping_status=UNRESOLVED_PUBLISHED_MEDIAN_RECONSTRUCTION`, but excluded from technical repeatability endpoints. This is a source-integrity exclusion, not an outcome-driven performance exclusion.

The validated set contains 10 sites across all 7 TFs. Every validated sequence has positive signal and at least 9 spots; TBP arrays contain 9 spots per sequence, and the remaining validated blocks contain at least 10 spots per sequence. No missing or non-positive spot signal was observed in the validated rows.

## Processing and split

The primary split is frozen by `V1_3_2_EXTERNAL_VALIDATION_CONTRACT.md`: sorted keys, one NumPy `PCG64` generator initialized at seed 1301, and alternating assignment within each `(TF, site, sequence, array)` stratum. The source workbook does not supply an array or batch identifier; processed rows therefore use `array_id=unknown_array` and record this limitation explicitly. The source parser retains raw signal and source identifiers. Split-specific effects are computed as log2(mean mutant signal / mean matched reference signal) within the split. No value is imputed.

The primary split and 100-seed secondary stability table are in `results/v1_3/external_samba_primary_split_effects.parquet` and `results/v1_3/external_samba_split_stability.tsv`.

## Missing or non-evaluable items

- GEO raw files are downloaded and retained, but the primary parsed endpoint uses the calibration workbook because it contains the documented site, position, pairing, and perturbation metadata needed for the canonical/mismatch key.
- Ets1 site 2 and site 4 are not used for endpoint summaries because their published medians cannot be reconstructed from the joined raw rows without an undocumented block-selection rule.
- These data are technical spot repeats. They cannot identify biological replicate variance or serve as an independent biological noise ceiling.
