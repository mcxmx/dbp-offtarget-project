# Competition assay semantics

The paper defines the plotted quantity as relative binding activity, PE/FITC normalized to the no-competitor condition. A nonbiotinylated competitor reduces the labeled-target binding signal. Therefore lower normalized PE/FITC means stronger competition by that mutant competitor.

The frozen experimental direction used in v1.1 is:

`experimental_effect = - normalized_PE/FITC`

Thus larger experimental effect means stronger competitor binding / stronger competition. This is a monotonic sign reversal only; no per-protein scaling, register selection, or model-driven normalization is applied. The raw normalized signal is preserved in the tidy table.

The XLS has heatmap means of two replicates, not separate replicate columns. Correlations therefore use the published two-replicate mean values and are labeled accordingly.
