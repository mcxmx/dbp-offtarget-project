# v1.3.2 External Data Provenance

Retrieval date: 2026-09-17. New source files were stored below `data/raw/v1_3_2_external/`. SHA256 values were computed once and persisted in `artifacts/cache/file_manifest.json`; subsequent use is size/mtime-aware.

| dataset | local source | official source | size | SHA256 | use |
|---|---|---|---:|---|---|
| Afek et al. Supplementary Table 1 (MOESM4) | `data/raw/v1_3_2_external/samba/41586_2020_2843_MOESM4_ESM.xlsx` | [Springer Nature source file](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-020-2843-2/MediaObjects/41586_2020_2843_MOESM4_ESM.xlsx) | 7,775,515 | `ba3afc5b182a319bd6b200b69b58996fde593af2712bb18393d4448951932c18` | SaMBA library/site inventory and raw mutation-array sheets |
| Afek et al. Supplementary Table 4 (MOESM7) | `data/raw/v1_3_2_external/samba/41586_2020_2843_MOESM7_ESM.xlsx` | [Springer Nature source file](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-020-2843-2/MediaObjects/41586_2020_2843_MOESM7_ESM.xlsx) | 12,968,460 | `09715f8531f2afffb81348e1d369e60b6d861fc4fb75e254693979ff2d1bd6a7` | 12 calibration sites, canonical/mismatch metadata, and raw spot-level signals |
| GEO GSE156375 raw archive | `data/raw/v1_3_2_external/samba/GSE156375_RAW.tar` | [NCBI GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156375) / [NCBI FTP archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156375/suppl/GSE156375_RAW.tar) | 15,114,240 | `a23e11e027ae34f1eed117bc8e3cad158085b0a7d761b32b046d9f4dfd166bf6` | Official raw sample archive retained for provenance; primary calibration parse uses MOESM7 metadata |

The primary SaMBA analysis is based on official publisher source data, not figure OCR. The parser records source rows, sequence, source sheet, spot IDs, raw signal, and the published median used for source-integrity validation. Array/batch identifiers are absent from the calibration workbook; processed rows therefore use `array_id=unknown_array` and explicitly record `array_id_source=not_supplied_by_source_workbook`.

## T227 provenance and recovery status

The primary official source is the SAMPDI-3D article and its PMC-hosted supplementary PDF: [Bioinformatics article DOI 10.1093/bioinformatics/btab567](https://doi.org/10.1093/bioinformatics/btab567), [PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10186157/). The supplement was inspected for row-level T227 data and per-row SAMPDI-3D/FoldX predictions; it provides only benchmark descriptions and aggregate values. No valid T227 row-level source file was downloaded into the repository. The former Clemson SAMPDI-3D endpoint was unavailable at audit time. See `reports/v1.3/V1_3_2_T227_DATA_AUDIT.md`.
