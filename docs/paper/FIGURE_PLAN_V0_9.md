# Figure Plan V0.9

## Figure 1: Benchmark and exposure-controlled design

Show the seven-protein GSE237017 cohort, protein-cluster leave-one-out split, canonical reverse-complement classes, development exposure boundary and protein-level statistical unit. Do not label this cohort external validation.

## Figure 2: Unseen-protein prediction performance

Show per-protein Spearman points for sequence k-mer, v0.5 M0-M3 and corrected v0.6 M0-M3, with median reference lines and seed stability. Do not use a median-only bar plot. The experimental replicate reference is a context line, not a ceiling.

## Figure 3: Conditioning collapse

Show protein-shuffle correlation, target-shuffle correlation and conditioning effect sizes for v0.5/v0.6 M3. Include the synthetic conditionality control as a positive control. Make the collapse threshold 0.995 explicit.

## Figure 4: Falsification experiments

Show the v0.5 RC collision audit and v0.6 injectivity repair, dense-supervision pilots and local/residue attention diagnostics. The visual message is that specific implementation, pair-density and local-attention explanations were tested and did not restore real-data conditionality.

## Figure 5: Hard cases and natural-designed shift

Show the 1,515 hard-case taxonomy, the 871/1,515 all-model-failure subset, and assay/construct/protein-distribution differences between natural PBM and designed uPBM. Mark causal interpretations as hypotheses. Structure-method coverage is annotated as 2/7 and not used as a universal failure panel.

## Figure 6: Independent external validation

Reserved. Do not draw a placeholder performance panel or reuse GSE237017. Once a Tier 1 or Tier 2 cohort passes the metadata firewall, show preregistered per-protein endpoints, coverage, failures and protein-level uncertainty. Tier 3 results belong in an orthogonal supplement unless the preregistered endpoint explicitly supports a primary panel.

## Supplementary structure diagnostic

Place DeepPBS and NA-MPNN per-protein points, their 2/7 coverage, PWM/PPM proxy definitions, and AF3/ContactSeek feasibility in supplementary material. They must not appear as a comprehensive SOTA ranking or as evidence that structure-based methods fail.
