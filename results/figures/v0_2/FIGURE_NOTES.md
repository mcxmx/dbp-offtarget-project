# Figure Notes for Benchmark v0.2

Retrieval/audit date: 2026-08-29

## fig1_dataset_overview_v0_2.png

Can show: the size of the curated v0.2 perturbation benchmark per PDB-derived target.
Cannot show: binding specificity or off-target risk.
Category: dataset result / pipeline validation.

## fig2_gc_distribution_v0_2.png

Can show: GC-matched controls track target GC content more closely than fully random DNA.
Cannot show: binding preference or protein-conditioned specificity.
Category: dataset result / pipeline validation.

## fig3_single_mutation_sequence_proxy_landscape_v0_2.png

Can show: how sequence-only proxy metrics change under single substitutions for one curated target.
Cannot show: biological specificity landscape or mutation effect on binding.
Category: sequence-only baseline.

## fig4_double_mutation_sequence_proxy_landscape_v0_2.png

Can show: how sequence-only proxy metrics change under double substitutions for one curated target.
Cannot show: epistasis, binding energy, or off-target risk.
Category: sequence-only baseline.

## fig5_sequence_proxy_distribution_v0_2.png

Can show: separation induced by sequence similarity between target, mutants, GC controls, and random controls.
Cannot show: protein-DNA binding specificity or calibrated safety margin.
Category: sequence-only baseline / pipeline sanity check.

## Shared Limitation

No v0.2 figure should be interpreted as a biological specificity result. The current scoring is not protein-conditioned and has no quantitative experimental specificity ground truth attached.
