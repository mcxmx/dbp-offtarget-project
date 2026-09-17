# Assay-Specific Validation Rules

These endpoints are frozen before quantitative labels are accessed. An assay is not forced into the GSE237017 7-mer Spearman estimand when its data structure does not support that interpretation.

## Tier 1: complete or near-complete landscape

Eligibility: complete 7-mer/8-mer landscape or a PBM-like assay with a defined register. Primary endpoint is median per-protein Spearman between predicted and measured full-landscape rankings. Secondary endpoints are per-protein top-k enrichment, NDCG at predeclared k, and hard-case recovery if a frozen hard-case set exists. Ambiguous register or unavailable construct sequence excludes a protein before labels are inspected.

## Tier 2: dense randomized library

Eligibility: many measured variants per protein with a documented randomized region and a deterministic mapping from reads to sequences. Primary endpoint is median per-protein rank correlation when the measured table supports quantitative ranks; otherwise use a predeclared discrimination metric (on-target/near-target enrichment or top-quantile recovery) selected from the schema before labels. Secondary endpoints are on-target versus near-target enrichment, top-quantile recovery, and motif-level agreement. The endpoint is selected by assay schema, never by result quality.

## Tier 3: sparse mutation or competition assay

Eligibility: finite tested variants or competition pairs. Primary endpoint is directional agreement of mutational effects or rank correlation over tested variants, as determined by the preregistered table schema. Secondary endpoint is single-base specificity recovery or pairwise competitor ordering. Tier 3 is orthogonal evidence only and cannot be reported as full-landscape validation.

## Shared rules

Protein is the independent unit. Report all eligible proteins, failures and coverage. Reverse-complement handling uses canonical RC classes. DNA register, strand convention, scoring rule, method set, checkpoint, and seed aggregation are fixed before labels. No post-hoc model, protein, metric or register selection is permitted.
