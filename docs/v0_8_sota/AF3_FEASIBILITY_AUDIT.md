# AF3 Feasibility Audit

A full 7-protein x 8,192-candidate complex landscape would require 57,344 separate protein-DNA complex predictions before contact extraction, with additional sensitivity to stochastic inference and DNA register. AF3 contact probability is a structural/contact diagnostic, not a calibrated affinity or PBM E-score predictor. ContactSeek's public workflow extracts contact probabilities from AF3 confidence JSON/embedding outputs and demonstrates an exposed Glasscock-designed example; it does not provide a frozen, independent seven-protein full landscape for this project.

No defensible cheap approximation was accepted: replacing per-candidate complex inference with a sequence-only surrogate would change the estimand and violate the structure-method audit. `af3_candidate_manifest.csv` therefore records a preregistered subset definition only (on-target, Hamming-1/2, the frozen hard-case set, and fixed-seed stratified background); `af3_contact_metrics.csv` is an explicit not-run schema. No AF3 result was used for method selection.

Status: `FEASIBILITY_ONLY`, not a performance benchmark and not external validation.
