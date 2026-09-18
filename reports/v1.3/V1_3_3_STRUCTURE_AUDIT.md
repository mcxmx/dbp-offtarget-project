# v1.3.3 Structure Audit

Status: frozen before natural-TF model/label performance inspection (2026-09-18).

The audit uses only RCSB structure metadata and SaMBA sequence/mutation-coordinate metadata. No SaMBA effect column is loaded.

## Selection rule

Primary candidates were fixed label-blind by TF identity and direct experimental DNA-bound annotation. The mapping uses the contract's deterministic exact/global alignment rule; no correlation, rank, or pairwise metric was inspected.

## Site audit

| TF | Site | PDB | Protein chain | DNA chains | Method | Resolution A | Motif | PDB DNA | Identity | Mapping | Reason |
|---|---|---|---|---|---|---:|---|---|---:|---|---|
| Cbf1 | Cbf1 | 8OVW | A | D,E | ELECTRON MICROSCOPY | 3.4 | TCGCACGTGCCA |  | NA | MAPPING_AMBIGUOUS | 2 tied label-blind alignments |
| Egr1 | Egr1 | 1AAY | A | B,C | X-RAY DIFFRACTION | 1.600 | GGCGTGGGCGAT |  | NA | MAPPING_AMBIGUOUS | no candidate alignment maps every mutated coordinate |
| Ets1 | Ets1_site_1 | 2STT | C | A,B | SOLUTION NMR |  | GATATCCGGTT | TCGAGCCGGAAGTTCGA | 0.727 | ELIGIBLE | unique label-blind mapping |
| Ets1 | Ets1_site_2 | 2STT | C | A,B | SOLUTION NMR |  | ATGTTCCGGTT | TCGAGCCGGAAGTTCGA | 0.727 | ELIGIBLE | unique label-blind mapping |
| Ets1 | Ets1_site_3 | 2STT | C | A,B | SOLUTION NMR |  | AGGGGTAAGCG | TCGAGCCGGAAGTTCGA | 0.727 | ELIGIBLE | unique label-blind mapping |
| Ets1 | Ets1_site_4 | 2STT | C | A,B | SOLUTION NMR |  | AAAAACCGCAG | TCGAGCCGGAAGTTCGA | 0.636 | ELIGIBLE | unique label-blind mapping |
| GR | GR | 1R4R | A | C,D | X-RAY DIFFRACTION | 3.0 | AGAACAGCATGTACA |  | NA | MAPPING_AMBIGUOUS | 2 tied label-blind alignments |
| Max | Max_site_1 | 1HLO | A | C,D | X-RAY DIFFRACTION | 2.800 | GATCACGTGAAT |  | NA | MAPPING_AMBIGUOUS | no candidate alignment maps every mutated coordinate |
| Max | Max_site_2 | 1HLO | A | C,D | X-RAY DIFFRACTION | 2.800 | GATCGCGTGAAT |  | NA | MAPPING_AMBIGUOUS | no candidate alignment maps every mutated coordinate |
| TBP | TBP_site_1 | 1TGH | A | B,C | X-RAY DIFFRACTION | 2.900 | CCTTTTATAG |  | NA | MAPPING_AMBIGUOUS | 4 tied label-blind alignments |
| TBP | TBP_site_2 | 1TGH | A | B,C | X-RAY DIFFRACTION | 2.900 | ACTTTTATAG |  | NA | MAPPING_AMBIGUOUS | 4 tied label-blind alignments |
| p53 | p53 | 1TUP | A | E,F | X-RAY DIFFRACTION | 2.200 | ACATGCCCGGGCATGC |  | NA | MAPPING_AMBIGUOUS | 2 tied label-blind alignments |

## Coverage

- SaMBA sites audited: 12 across 7 TFs.
- Eligible mappings: 4 sites.
- Excluded/ambiguous mappings: 8 sites.
- Primary v1.3.2 source-integrity-valid sites among the structure-eligible set: 2 (Ets1_site_1 and Ets1_site_3). Ets1_site_2 and Ets1_site_4 are retained in the raw audit but excluded from primary labels because their published medians failed the frozen source reconstruction.
- `protein_sequence_coverage` is reported conservatively as structure-chain coverage metadata; no sequence identity to the assay construct was inferred without a frozen reference mapping.
- If the eligible count is too small for a TF-level challenge, the result is reported as `STRUCTURAL_COVERAGE_LIMITED` and v1.4 remains closed.
