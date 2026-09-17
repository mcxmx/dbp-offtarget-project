# Competition source-data audit

Source: `data/raw/v1_1_competition/source_data_extended_data_fig3.xls`

- Public PMC URL: https://pmc.ncbi.nlm.nih.gov/articles/instance/12618268/bin/41594_2025_1669_MOESM16_ESM.xls
- Resolved download URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM16_ESM.xls
- DOI: 10.1038/s41594-025-01669-4
- Retrieval date: 2026-09-14
- File size: 72192 bytes
- SHA256: `5ba651f551af207beb1fb91be255317242bd93b213d5b35de7a032cf1b76ace0`
- Workbook engine: xlrd/OLE XLS
- Sheets: 10
- Relevant original designs parsed: DBP001, DBP003, DBP005, DBP006, DBP009, DBP035, DBP048

Each sheet has a five-column table (`position`, `original_base`, `new_base`, `sample`, and `Median PE/FITC (Normalized)`). Relevant sheets contain 13 or 14 positions with three non-WT substitutions per position. DBP005, DBP009, and DBP035 include an additional WT row; the other relevant sheets do not.

No merged cells or missing values were observed. The sheets contain mutation rows (and WT rows in three relevant sheets) but no explicit no-competitor control row; the provided signal is already normalized to the no-competitor condition. There are no replicate columns. The XLS values are the heatmap/source values described by the paper as the mean of two replicates. Individual replicate values are not present in this downloaded workbook, so the tidy records are labeled `mean_of_two_replicates`; no replicate-level values are invented.

The source includes DBP023, DBP056, and DBP062 as additional designs. They are retained in the inventory audit but excluded from the seven-protein analysis.
