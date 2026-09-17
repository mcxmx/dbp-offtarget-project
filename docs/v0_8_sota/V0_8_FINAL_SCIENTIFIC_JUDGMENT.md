# V0.8 Final Scientific Judgment

1. **Is the gap only in sequence models?** Current evidence supports a sequence-model gap. It does not test enough independent structure predictions to claim the gap is universal across computational methods.
2. **Do DeepPBS/NA-MPNN shrink it?** No demonstrated rescue: DeepPBS median 0.1593 and NA-MPNN -0.0407, each on only two exposed proteins. Unequal coverage and N=2 make this descriptive, not a population comparison.
3. **Which method is closest to reproducibility?** The experiment-vs-experiment reference (median 0.5914) is the closest reference by definition; no model is close on a valid seven-protein, matched protocol. Among evaluated models, the sequence k-mer baseline (0.2321) exceeds the partial structure medians, but this is not evidence that it is biologically superior.
4. **Are 871/1515 common hard cases also structure hard cases?** Cannot be answered. Only 398 hard candidates in DBP35/DBP48 are structure-evaluable; the remainder have no legal structure prediction.
5. **Title wording:** use **sequence-based specificity prediction**. "Computational specificity prediction" would overclaim until a sufficiently covered structure benchmark is completed.
6. **Independent cohort:** none is currently qualified. RFdiffusion3 is an unexposed metadata candidate, not a validation cohort; labels remain locked.
7. **Pre-submission blocker:** `BLOCKED ON INDEPENDENT CONFIRMATORY DATA` and, secondarily, structure coverage for a fair SOTA comparison.
8. **Q1 benchmark paper evidence:** sufficient for a carefully scoped Route B benchmark/generalization-gap manuscript, not for a universal mechanism claim or Route A predictive-method claim.

**Strongest claim:** Under strict unseen-protein evaluation, existing sequence-based approaches show a reproducible benchmark-level generalization gap on the present de novo DBP cohort; representation repair improves ranking but does not restore protein- or target-dependent prediction.
