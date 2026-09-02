# Natural PBM Control Plan for v0.5+

The v0.4 designed-DBP benchmark uses processed uPBM E-scores from GSE237017. A future natural-to-designed OOD claim needs assay-matched controls because a natural HT-SELEX to designed uPBM comparison would confound protein distribution shift with assay distribution shift.

Priority sources to investigate next:

- CIS-BP PBM-derived specificity tables where protein IDs and protein sequences are resolvable.
- UniPROBE / universal PBM datasets with raw or processed k-mer enrichment and protein identity.
- JASPAR entries backed by PBM/uPBM experiments, used only when the experimental provenance is clear.

Recommended comparisons:

- Natural PBM train to natural PBM held-out protein.
- Natural PBM train to designed GSE237017 uPBM external test.
- Protein-family-out splits within natural PBM.
- Designed leave-one-protein/cluster-out only as secondary, because N=7 is small.

All DNA splits must group reverse-complement equivalent sequences together.
