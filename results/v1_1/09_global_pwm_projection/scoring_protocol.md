# v1.1 global projection protocol

Candidate universe: 8,192 canonical reverse-complement 7-mers. For DeepPBS, the fixed mapped seven model positions are scored as `sum(log(P+1e-12))`; forward and reverse-complement candidate orientations are both scored and the maximum is retained. NA-MPNN uses its frozen v1.0 canonical landscape verbatim because native PPM tensors were not retained. No positionwise RC averaging and no PBM-driven register optimization are used.
