# Manuscript Outline

## Introduction

Aggregate DNA-binding specificity metrics do not distinguish protein conditionality, position sensitivity, and within-position nucleotide identity. The benchmark decomposes a mutation effect into a position component and a nucleotide-specific residual while treating assay estimands as distinct observables.

## Results

### R1. Protein-conditioned global models exhibit conditionality collapse

Retain the frozen v0.5-v0.6 protein/target shuffle diagnostics as an orthogonal protein-side axis. These are development-exposed diagnostics, not independent validation.

### R2. Designed-DBP mutation prediction is dominated by position-level information

The frozen v1.3 development benchmark separates global, position, and identity residual performance. Overall and position-level agreement are materially higher than within-position identity agreement for DeepPBS.

### R3. Prospective locked holdout reproduces position >> identity

DBP023/056/062 were predicted before label access and then evaluated under the frozen contract. All three directionally reproduce position-dominant performance; the holdout n is three.

### R4. Independent SaMBA data quantify technical identifiability

Canonical Watson-Crick SaMBA spot-level technical splits show high global, position, and within-position identity repeatability. This is technical within-assay repeatability, not biological reproducibility or a cross-assay noise ceiling. Mismatches are reported separately.

### R5. Natural-TF mutation benchmark status

The official SAMPDI-3D T227 supplement exposes only aggregate benchmark values and contains no auditable row-level table or predictions. The natural-TF decomposition is therefore not evaluable in v1.3.2, and no generalization claim is made.

## Discussion

- Designed-protein n remains small and development-exposed except for the three locked holdouts.
- Technical spot repeats must not be described as biological replicates.
- PBM, SaMBA, and competition assay values have distinct estimands and context.
- The external SaMBA result argues against treating identity information as universally unobservable, but does not remove assay-specific uncertainty in competition data.
- DeepPBS-specific identity failure remains supported by the frozen designed benchmark but awaits an independent natural-TF row-level benchmark for broader generalization.
- Future benchmarks should require replicate-resolved mutation matrices and pre-registered global/position/identity decomposition.

## Current model-development decision

Do not open v1.4 model development until a row-level natural-TF benchmark is recovered or a new independent assay is established under the frozen contract. No GNN is trained in v1.3.2.
