# NA-MPNN representation note

The frozen v1.0 artifact contains a canonical 7-mer score landscape, not the native `predicted_ppm` tensor. The conceptual PWM table is therefore emitted with NA probabilities; local NA-MPNN mutation effects are labeled as derived from the frozen landscape, never as native PPM values.
