# v1.3.4 Technical Rescue Protocol

Status: frozen before any post-reveal rescue attempt.

## Scope

This protocol preserves the v1.3.3 natural-TF primary result. The frozen
challenge remains PDB 2STT, ETS1 chain C with DNA chains A/B. Any later
prediction is labelled POST_REVEAL_TECHNICAL_SENSITIVITY; it cannot become
prospective validation, primary validation, or an independent blind benchmark.
No new GNN, neural network, fine-tuned model, checkpoint, or training run is
allowed in v1.3.4.

## Deterministic structure rule

2STT is a multi-model NMR ensemble. The primary rescue structure is
data/processed/v1_3_4/2STT_model1.pdb, extracted from MODEL 1 only. The
operation is deterministic extraction plus necessary PDB formatting only. It
does not perform minimization, geometry optimization, DNA remodeling, mutation,
coordinate fitting, or label-guided register adjustment. Protein chain C and
both complementary DNA strands are retained with residue numbering and
coordinates. MODEL 13 and MODEL 25 may be generated only for predeclared
prediction-vs-prediction stability, never selected by experimental performance.

## Bounded rescue routes

1. Positive control and Route 1: verify the existing successful DBP035 frozen
   DeepPBS artifact and reuse the exact local/frozen pipeline, environment,
   commit, checkpoints, and preprocessing for 2STT_model1.pdb.
2. Route 2: only if Route 1 fails, attempt the official documented Docker
   pipeline. Record DeepPBS version/commit, checkpoint, image, command,
   environment, hardware, and CPU/GPU mode. A non-identical implementation is
   IMPLEMENTATION_SENSITIVITY and cannot replace the frozen primary result.
3. Route 3: the official webserver is sanity-check only. Without exact
   model, checkpoint, preprocessing, and version provenance it is
   WEB_SERVER_DESCRIPTIVE_ONLY and excluded from quantitative comparison.

The rescue has exactly three bounded actions: positive control; local/frozen
MODEL 1; official Docker MODEL 1 if local fails. Do not try alternate PDBs,
registers, orientations, signs, checkpoints, preprocessing variants, or
performance-guided debugging. Stop after Route 1/2 regardless of outcome.

## Output and interpretation

If MODEL 1 produces an NPZ, write
results/v1_3/natural_tf_2stt_model1_posthoc_predictions.tsv and the frozen
global, position, identity-residual, and pairwise metrics to
natural_tf_2stt_model1_posthoc_decomposition.tsv. These outputs must state
POST-REVEAL TECHNICAL SENSITIVITY, one TF, limited site count, and that labels
were revealed before rescue. They are not prospective validation. If no NPZ is
produced, retain the existing NOT_EVALUABLE primary result and record the
actual failure stage and exit code.

## Reproducibility and caching

The persistent artifacts/cache/file_manifest.json is authoritative. Hashing
is size/mtime-aware and is performed once for a newly created artifact only.
Existing files are never rehashed when size and mtime are unchanged. Existing
prediction tables, SaMBA parses, and figure source tables are reused directly
when their input/config/model/mapping provenance is unchanged.
