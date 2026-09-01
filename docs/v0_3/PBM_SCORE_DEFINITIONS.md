# PBM Score Definitions for GSE237017

Source: GEO accession GSE237017, sample-level processed uPBM files, and GEO data processing text.

## Experiment

GSE237017 contains universal protein-binding microarray experiments for designed DNA-binding proteins DBP1, DBP3, DBP5, DBP6, DBP9, DBP35, and DBP48. GEO describes the arrays as 15K Agilent dsDNA arrays designed so that all possible 9-bp sequences are covered, with every 7-mer represented in at least 16 spots.

## Data Processing Recorded by GEO

The GEO sample metadata states that Alexa 488 signals were position-adjusted to correct microarray non-uniformity. A Seed-and-wobble algorithm was then used to compute enrichment scores and median intensities for all possible 7-mers. GEO defines the processed values as E-scores, median intensities, and z-scores for 7-mers.

## E-score

Recommended v0.3 primary score: `PBM E-score`, stored as `experimental_score_primary`.

Interpretation:

- Experimental PBM specificity/enrichment score for 7-mer ranking within a protein.
- Useful for per-protein specificity landscape analysis.
- Suitable as the first primary target for designed-DBP external benchmark ranking.

Not justified:

- Not Kd.
- Not binding free energy.
- Not binding probability.
- Not directly comparable as absolute affinity across different proteins.

## Median Intensity

Interpretation:

- Processed microarray signal intensity summary for 7-mers.
- Useful as a secondary assay-derived measurement.

Not justified:

- Not directly calibrated affinity.
- Not automatically comparable across proteins or experiments without normalization.

## Z-score

Interpretation:

- Standardized processed score included by GEO for 7-mers.
- Useful as a secondary consistency/QC measure.

Not justified:

- Not a probability or Kd.
- Not used as the primary v0.3 score unless future analysis shows it is preferable.

## Reverse Complement Handling

The processed 7-mer files contain two 7-mer columns. In v0.3 the parser explicitly expands both the primary and reverse-complement companion columns, assigning the same GEO-provided scores to both 7-mer orientations. `reverse_complement_qc.csv` confirms zero difference between paired reverse-complement scores after expansion.
