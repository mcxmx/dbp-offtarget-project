# Paper uPBM Reproduction Notes

Audit date: 2026-09-01

Primary source: Nature article "Computational design of sequence-specific DNA-binding proteins", DOI 10.1038/s41594-025-01669-4.

Additional sources:

- GEO GSE237017 sample metadata and processed uPBM files.
- Source Data Extended Data Fig. 8 (`41594_2025_1669_MOESM20_ESM.xls`).
- Source Data Fig. 4 (`41594_2025_1669_MOESM12_ESM.xlsx`).
- Supplementary Tables 1-3 (`41594_2025_1669_MOESM3_ESM.xlsx`).

## Paper Definition Captured

The paper reports that uPBM E-scores were used to evaluate whether 7-mers containing the designed binding-site motif were enriched among high-scoring sequences. Extended Data Fig. 8 reports motif percentile values for DBP6, DBP9, DBP48, DBP5, DBP35, DBP1, and DBP3.

The Methods/GEO processing notes describe processed PBM E-scores computed with Seed-and-wobble from Alexa 488 signal after position adjustment. These are processed experimental uPBM specificity/enrichment scores, not Kd, free energy, or binding probability.

## v0.3.1 Reproduction Rule

The Extended Data Fig. 8 source-data workbook stores 8192 rows. Each row has two reverse-complement 7-mer columns and one E-score per replicate. v0.3.1 therefore treats each source-data row as one reverse-complement equivalence class.

For each DBP:

1. Rank the 8192 source-data rows within each replicate by E-score.
2. Convert rank to percentile, with higher E-score giving higher percentile.
3. Select rows where either 7-mer column contains the designed motif.
4. Compute the mean percentile of selected rows per replicate.
5. Average replicate-level motif percentiles when replicates are available.

Tolerance: <= 2 percentile points absolute difference from the paper-reported value.

## DBP48 Target Clarification

Supplementary Table 1 records DBP48 as originally designed against target `I_b`, whose top strand is `CGCCCAAAGCCGCG`. The Fig. 4 caption states that DBP48 was analyzed with sequence C because of improved binding signal and nearly identical modeled binding sites. v0.3.1 therefore records:

- Original design target: target I.
- Experimental assay/PBM evaluation reference: target C.
- PBM motif used for percentile reproduction: `CTGACG`.

This separates original design target, assay target, and PBM motif reference.

## Reproduction Result

Overall status: PASS

| protein_id | paper_reported_percentile | our_reproduced_percentile | absolute_difference | motif_sequence | motif_length | n_matching_7mers | reproduction_status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DBP1 | 33.190 | 33.576 | 0.386 | GCAGG | 5 | 67;67 | PASS |
| DBP3 | 46.580 | 47.834 | 1.254 | GCAGGA | 6 | 8;8 | PASS |
| DBP35 | 81.880 | 82.615 | 0.735 | TGCACA | 6 | 8;8 | PASS |
| DBP48 | 97.590 | 97.950 | 0.360 | CTGACG | 6 | 16 | PASS |
| DBP5 | 86.540 | 86.026 | 0.514 | TGCACA | 6 | 8 | PASS |
| DBP6 | 99.540 | 99.252 | 0.288 | TGCACA | 6 | 8;8 | PASS |
| DBP9 | 99.890 | 99.807 | 0.083 | TGCACA | 6 | 8;8 | PASS |

## Direct Paper-Derived Versus Implementation Choices

Direct from paper/source data:

- GSE237017 is the uPBM accession for the designed DBPs.
- Extended Data Fig. 8 reports the DBP-specific motif percentiles.
- The source-data workbook provides replicate E-scores for 8192 7-mer reverse-complement rows.
- Supplementary Table 1/2 provides design target IDs and exact dsDNA target sequences.
- Fig. 4 describes DBP48 analysis with sequence C.

Implementation choices:

- Use the source-data row as the RC-equivalence class.
- Match motif if either 7-mer column in that row contains the motif.
- Average replicate-level motif percentiles for DBPs with two replicates.
- Use the motif mapping recorded in `metadata/v0_3_1/designed_dbp_target_definitions.csv`.
