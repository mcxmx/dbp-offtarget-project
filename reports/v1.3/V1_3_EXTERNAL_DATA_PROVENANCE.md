# v1.3.1 External Data Provenance

Retrieval date: 2026-09-17. Paper DOI: [10.1038/s41594-025-01669-4](https://doi.org/10.1038/s41594-025-01669-4). Files were downloaded from the publisher's official Springer Nature static-content host. SHA256 was calculated once after download and cached in `artifacts/cache/file_manifest.json`; subsequent checks use path, size, and mtime first.

| figure | original filename | publisher source | bytes | SHA256 | contents used | replicate values | transformation |
|---|---|---|---:|---|---|---|---|
| Extended Data Fig. 3 | `41594_2025_1669_MOESM16_ESM.xls` | [official workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM16_ESM.xls) | 72192 | `5ba651f551af207beb1fb91be255317242bd93b213d5b35de7a032cf1b76ace0` | 10 DBP competition landscapes | No: one published summary column only | Mutation rows parsed; `normalized_effect = -raw_measurement`; prior mean parity checked |
| Extended Data Fig. 9 | `41594_2025_1669_MOESM21_ESM.xls` | [official workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM21_ESM.xls) | 15216640 | `8be665a7c8cc577e93f883407e5cdf478ab4a87c809a21e7612a033c2092a7db` | panel 9e DBP35 K18V/R33N/P42Q landscape | No replicate-valued columns in panel 9e | Mutation rows parsed; `normalized_effect = -raw_measurement`; other raw-event sheets retained untouched |

The Extended Data Fig. 3 sheet names unexpectedly use `Extended_Data_Figure_1_*`, but their ten DBP identifiers and values match the previously archived Extended Data Fig. 3 workbook. No screenshot digitization or OCR was used. The large Fig. 9 workbook is retained verbatim; only panel 9e was converted for this benchmark.
