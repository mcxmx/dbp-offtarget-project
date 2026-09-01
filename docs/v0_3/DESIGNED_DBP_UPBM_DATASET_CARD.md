# Designed DBP uPBM Dataset Card v0.3

## Dataset Name

Designed DBP uPBM experimental specificity benchmark v0.3.

## Source

- GEO accession: GSE237017
- GEO title: Computational design of sequence-specific DNA-binding proteins
- Paper reference: DOI 10.1038/s41594-025-01669-4; PMID 40940539
- Primary GEO URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237017

## Proteins

The benchmark contains 7 designed DNA-binding proteins: DBP1, DBP3, DBP5, DBP6, DBP9, DBP35, and DBP48.

Protein sequences recovered from official supplementary material: 7/7.
Intended target DNA sequences recovered from official supplementary material: 7/7.

## Samples and Replicates

Usable GEO samples: 12.
Replicate pairs with E-score QC: 5.
Proteins with replicate samples: DBP1; DBP3; DBP35; DBP6; DBP9.
Single-replicate-only proteins: DBP48; DBP5.

| gsm_id | protein_id | protein_concentration | replicate | sample_title |
| --- | --- | --- | --- | --- |
| GSM7593072 | DBP1 | 300uM | 1 | DBP1 at 300uM concentration replicate 1 |
| GSM7593079 | DBP1 | 300uM | 2 | DBP1 at 300uM concentration replicate 2 |
| GSM7593073 | DBP3 | 75uM | 1 | DBP3 at 75uM concentration replicate 1 |
| GSM7593080 | DBP3 | 75uM | 2 | DBP3 at 75uM concentration replicate 2 |
| GSM7593077 | DBP35 | 100uM | 1 | DBP35 at 100uM concentration replicate 1 |
| GSM7593083 | DBP35 | 100uM | 2 | DBP35 at 100uM concentration replicate 2 |
| GSM7593078 | DBP48 | 75uM | 1 | DBP48 at 75uM concentration replicate 1 |
| GSM7593074 | DBP5 | 300uM | 1 | DBP5 at 300uM concentration replicate 1 |
| GSM7593075 | DBP6 | 200uM | 1 | DBP6 at 200uM concentration replicate 1 |
| GSM7593081 | DBP6 | 200uM | 2 | DBP6 at 200uM concentration replicate 2 |
| GSM7593076 | DBP9 | 75uM | 1 | DBP9 at 75uM concentration replicate 1 |
| GSM7593082 | DBP9 | 75uM | 2 | DBP9 at 75uM concentration replicate 2 |

## Experimental Protocol

Universal protein-binding microarray experiments were performed for small recombinant designed DNA-binding proteins. The processed files provide 7-mer E-scores, median intensities, and z-scores.

## DNA Sequence Space

Each sample covers 16384 unique 7-mers after explicit expansion of the reverse-complement companion column. The final consensus benchmark contains 114688 protein-7mer measurements.

## Score Definitions

Primary score: PBM E-score, stored as `experimental_score_primary`.

Secondary scores: `median_intensity_mean`, `z_score_mean`.

The primary score is an experimental PBM specificity/enrichment score for per-protein ranking. It is not a Kd, binding probability, binding free energy, or absolute cross-protein affinity.

## Processing Steps

1. Downloaded GSE237017 family SOFT metadata and GSM supplementary files.
2. Parsed sample title, DBP ID, concentration, replicate, platform, and supplementary URLs from GEO metadata.
3. Downloaded all 12 processed 7-mer files and 12 raw spot-data files.
4. Preserved SHA256 and file size for every raw file.
5. Parsed processed 7-mer tables and expanded primary/reverse-complement 7-mer columns.
6. Built replicate-level consensus by mean E-score within each protein and 7-mer.
7. Added per-protein rank, percentile, and within-protein z-score without overwriting raw E-score.

## QC

- E-score replicate Pearson range: 0.551 to 0.765; median 0.665.
- E-score replicate Spearman range: 0.469 to 0.734; median 0.591.
- Reverse-complement score differences are zero after expansion.
- Sample coverage is complete for all processed samples after RC-column expansion.

## Intended Targets

| protein_id | intended_target_dna | target_length | target_id | target_context | confidence |
| --- | --- | --- | --- | --- | --- |
| DBP1 | TAGCAGGATGTGT | 13 | A | c | high |
| DBP3 | TAGCAGGATGTGT | 13 | A | c | high |
| DBP35 | GCAGATCTGCACATC | 15 | B | c | high |
| DBP48 | CGCCCAAAGCCGCG | 14 | I | b | high |
| DBP5 | GCAGATCTGCACATC | 15 | B | c | high |
| DBP6 | GCAGATCTGCACATC | 15 | B | b | high |
| DBP9 | GCAGATCTGCACATC | 15 | B | b | high |

## Known Limitations

- Processed specificity is 7-mer based; full target DNA binding is not directly measured by these tables.
- DBP5 and DBP48 have one GEO sample each in the parsed metadata, so replicate QC cannot be computed for them.
- Absolute E-scores should not be compared across proteins as binding affinity.
- The benchmark is in vitro uPBM, not in vivo genomic binding.

## Recommended Use

- Designed DBP specificity ranking.
- External validation for protein-conditioned DNA specificity models.
- Out-of-distribution evaluation of models trained on natural DBPs.

## Prohibited Interpretation

- Absolute affinity regression across proteins.
- Direct Kd or free-energy interpretation.
- Direct in vivo genomic binding claim.
- Claiming full-target affinity from overlapping 7-mer summaries.
