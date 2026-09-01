# v0.3.1 Validation Report

Audit date: 2026-09-01

Decision: GO TO V0.4

## 1. Paper PBM Motif Percentile Reproduction

Status: PASS

All 7 published Extended Data Fig. 8 motif percentiles were reproduced within the predeclared tolerance of 2 percentile points.

Maximum absolute difference: 1.2535 percentile points.

| protein_id | paper_reported_percentile | our_reproduced_percentile | absolute_difference | motif_sequence | motif_length | reproduction_status |
| --- | --- | --- | --- | --- | --- | --- |
| DBP1 | 33.190 | 33.576 | 0.386 | GCAGG | 5 | PASS |
| DBP3 | 46.580 | 47.834 | 1.254 | GCAGGA | 6 | PASS |
| DBP35 | 81.880 | 82.615 | 0.735 | TGCACA | 6 | PASS |
| DBP48 | 97.590 | 97.950 | 0.360 | CTGACG | 6 | PASS |
| DBP5 | 86.540 | 86.026 | 0.514 | TGCACA | 6 | PASS |
| DBP6 | 99.540 | 99.252 | 0.288 | TGCACA | 6 | PASS |
| DBP9 | 99.890 | 99.807 | 0.083 | TGCACA | 6 | PASS |

## 2. DBP48 Target Definition

DBP48 is separated into three concepts:

- Original design target: I_b / `CGCCCAAAGCCGCG`
- Experimental assay target: C / `CGACACCTGACGCG`
- PBM evaluation motif: `CTGACG`

This fixes the v0.3 ambiguity where a single `intended_target_dna` field could mix original design and assay/PBM evaluation references.

## 3. Reverse-Complement Units

Status: PASS

- Total oriented rows: 114688
- Oriented 7-mers per protein: 16384
- Total protein-RC-class units: 57344
- RC classes per protein: 8192

The confirmed independent sequence unit for v0.3.1 is the protein-RC-class unit, not the oriented row.

## 4. RC-Aware Sequence Baseline

The v0.3.1 baseline uses paper motifs and reverse-complement-aware comparison on 8192 RC classes per protein.

RC-aware hamming Spearman range: -0.194 to 0.284.
RC-aware 3-mer Spearman range: -0.139 to 0.362.

Median old-vs-new comparison:

| metric | old_spearman | rc_aware_spearman | delta |
| --- | --- | --- | --- |
| edit | 0.027 | 0.061 | -0.008 |
| hamming | 0.018 | 0.056 | 0.003 |
| kmer3 | 0.036 | 0.232 | 0.196 |
| kmer4 | 0.045 | 0.154 | 0.146 |

Interpretation: sequence-only similarity remains limited. Hamming/edit similarity is weak. Motif-level k-mer overlap explains part of the uPBM landscape for some DBPs, but it is still far below replicate agreement and is not protein-conditioned.

## 5. Disagreement Candidates

Status: PASS

The previous number 140 was the size of a per-protein top-20 examples table. It was not the total count.

v0.3.1 total sequence-vs-experiment disagreement candidates: 1515

Top examples table size: 140

Criterion: per-protein processed uPBM E-score >= 95th percentile and RC-aware motif hamming similarity <= protein median.

## 6. Replicate Noise Ceiling

This is an empirical replicate agreement / assay reproducibility reference, not a strict mathematical maximum.

| protein_id | pearson_correlation | spearman_correlation | n_aligned_7mers |
| --- | --- | --- | --- |
| DBP1 | 0.558 | 0.508 | 16384 |
| DBP3 | 0.739 | 0.655 | 16384 |
| DBP35 | 0.665 | 0.591 | 16384 |
| DBP6 | 0.551 | 0.469 | 16384 |
| DBP9 | 0.765 | 0.734 | 16384 |

E-score Pearson median: 0.665.
E-score Spearman median: 0.591.

## 7. Protein and Target Independence

Protein sequence clusters at 0.60 identity threshold: 4

Original target groups: 4
Assay target groups: 4
Motif groups: 4

Future splits must respect protein sequence clusters, target groups, motif groups, and canonical RC DNA classes.

## 8. DeepPBS / NA-MPNN Readiness

The dataset is ready to enter v0.4 as a benchmark arena for protein-conditioned baselines, including later DeepPBS/NA-MPNN-style comparisons if those tools are used only as evaluated scoring backends.

It is not ready for calibrated off-target risk claims, full-target affinity claims, or uncertainty calibration.

## 9. Remaining Limitations

- GSE237017 v0.3.1 uses processed uPBM 7-mer E-scores, not raw array-level reprocessing.
- Full designed targets are longer than 7 bp; motif percentile reproduction does not directly measure full-target affinity.
- DBP5 and DBP48 have single replicate only in the parsed GEO metadata.
- Natural-to-designed external evaluation can be confounded by assay shift unless natural PBM/uPBM controls are added.
- Sequence-only baseline is a sanity check, not a protein-conditioned model.

## 10. Gate Checks

| Requirement | Status |
| --- | --- |
| PBM parsing no major issue | PASS |
| Published motif percentile reproduction | PASS |
| RC handling explicit and tested | PASS |
| Target definitions separated | PASS |
| Disagreement count corrected | PASS |
| Tests all pass | PASS: 29 passed |

Final gate: GO TO V0.4
