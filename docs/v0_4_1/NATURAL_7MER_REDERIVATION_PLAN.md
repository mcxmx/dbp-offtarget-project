# Natural 7-mer Rederivation Plan

v0.4.1 does not rederive natural 7-mer scores from UniPROBE probe-level data.

Scientific reason: direct 8-mer to 7-mer truncation would merge multiple 8-mer contexts and change the PBM score definition. A valid 7-mer rederivation would need normalized probe-level data, probe sequence design, background model choices, and a documented PBM enrichment pipeline.

Future feasible plan:

1. Select one UniPROBE publication with downloadable deBruijn probe sequences and normalized probe intensity files.
2. Recompute 7-mer probe occurrence features with reverse-complement equivalence.
3. Define enrichment using a fixed robust statistic before looking at designed benchmark performance.
4. Compare rederived 8-mer scores to official UniPROBE 8-mer scores as a sanity check.
5. Only then export rederived 7-mer profiles.
