# External Cohort Selection Protocol

This protocol is written before quantitative label access. A confirmatory cohort must contain de novo designed DNA-binding proteins, construct/protein sequences, intended targets, quantitative specificity measurements, multiple independently designed proteins (preferably >7), and an assay compatible with the current canonical 7-mer ranking task. Complete 7-mer landscapes or a sufficiently dense randomized library are preferred.

Metadata inspection may establish study existence, assay type, sequence availability and raw-data availability. It may not inspect per-sequence scores, performance labels, easy/hard protein annotations, or any result used for model selection. GSE237017 and ContactSeek's Glasscock example are exposed and excluded from confirmation.

A cohort is not "qualified" until all inclusion fields are verified without labels, the frozen validation plan is committed, and the manifest is hashed. If no cohort meets the criteria, report `NO SUITABLE PUBLIC CONFIRMATORY COHORT IDENTIFIED` and request data from authors.
