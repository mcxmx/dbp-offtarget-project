# GO / NO-GO for v0.3 Designed DBP uPBM Benchmark

Decision: CONDITIONAL GO

## 1. Is GSE237017 sufficient as a designed-DBP external benchmark?

Yes, conditionally. It provides 12 usable uPBM samples across 7 designed DBPs and 114688 protein-7mer measurements with GEO provenance. It is suitable for external/OOD ranking evaluation, not absolute affinity regression.

## 2. Are replicates sufficiently consistent?

Partially. E-score replicate Pearson median is 0.665 and Spearman median is 0.591. This is usable for ranking analyses with QC caveats. DBP5 and DBP48 have single replicate only in GEO metadata.

## 3. Were all designed DBP protein sequences recovered?

Yes. 7/7 protein sequences were recovered from the official Nature supplementary workbook.

## 4. Were intended targets recovered?

Yes. 7/7 intended target DNA sequences were recovered from the official Nature supplementary workbook.

## 5. How does sequence-only similarity correlate with experimental specificity?

Weakly. Per-protein sequence-only Spearman correlations range from -0.049 to 0.221, with median 0.040. This supports using v0.3 to test protein-conditioned models.

## 6. Are there high-score 7-mers not explained by simple sequence similarity?

Yes. The v0.3 analysis found 140 candidate rows with top 1% PBM E-score and hamming similarity to target-derived 7-mers at or below the protein median.

## 7. Is the project ready for protein-conditioned baselines?

Conditionally yes. The dataset is ready for per-protein ranking baselines and zero-shot designed-DBP evaluation. It is not ready for calibrated affinity or uncertainty claims.

## Main Limitation

The current benchmark is 7-mer in vitro uPBM specificity data. Full intended target sequences are longer than 7 bp, so target-rank summaries use overlapping 7-mers and must not be treated as full-target binding affinity.

## Recommended Next Baseline

Run a protein-conditioned but non-neural baseline first: encode the designed DBP sequence and candidate 7-mer, evaluate per-protein Spearman and top-k enrichment, and compare against the v0.3 sequence-only baseline.
