# Model Split Plan for v0.4+

GSE237017 designed DBPs should be protected as an external designed-protein benchmark rather than mixed into natural transcription-factor training. v0.3.1 further requires reverse-complement equivalence and protein/target grouping controls.

## Split A: Natural Protein Random Split

Use only natural DBP experimental specificity datasets. Randomly split proteins after deduplicating identical protein sequences and highly similar motifs. This split is useful for debugging but should not be the headline generalization result.

## Split B: Protein-Family Split

Cluster natural proteins by sequence identity, domain annotation, and family labels where available. Hold out entire protein families to measure family-level generalization.

## Split C: Designed DBP Zero-Shot External Test

Train on natural DBPs only. Evaluate DBP1, DBP3, DBP5, DBP6, DBP9, DBP35, and DBP48 from GSE237017 as external/OOD designed DBPs. This is one of the key paper-level experiments: testing whether natural-trained protein-DNA models suffer an OOD performance drop on de novo designed binders.

Important caveat: if natural training data use HT-SELEX and designed external data use uPBM, observed performance drop may combine protein distribution shift and assay distribution shift. Assay-matched natural PBM/uPBM controls are needed before making a strong OOD claim.

## Split D: Leave-One-Designed-DBP-Out

After a zero-shot evaluation is established, optionally train/calibrate on six designed DBPs and test on the held-out designed DBP. This evaluates transfer within the designed-binder distribution, but it should not replace Split C.

## Leakage Controls

- Do not allow identical protein sequences across train/test.
- Cluster related natural proteins by family.
- Track DNA target similarity and motif similarity separately.
- Keep GSE237017 processed uPBM 7-mer E-scores separate from synthetic mutation benchmarks.
- Group train/test DNA by canonical reverse-complement equivalence class.
- Never place a sequence in train and its reverse complement in test.
- Do not treat DBP1/DBP3, DBP5/DBP35, or DBP6/DBP9 as fully independent strongest generalization pairs without cluster-aware analysis.
- Report per-protein metrics rather than pooling absolute scores across proteins.

## Recommended Metrics

- Per-protein Spearman correlation.
- Top-k enrichment within protein.
- Ranking quality for high-scoring 7-mers.
- OOD drop from natural-family held-out sets to designed DBPs.
