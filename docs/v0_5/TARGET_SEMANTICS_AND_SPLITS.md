# v0.5 Target Semantics and Split Contract

The v0.5 contract is generated from existing curated metadata by
`analysis/v0_5/00_build_target_manifest_and_splits.py`.

## Target semantics

`primary_target` is an alias for the independently reported
`intended_design_target`. It is sourced from the v0.3.1 target definitions and
is not inferred from the PBM landscape.

The manifest separately stores:

- `intended_design_target`
- `experimental_assay_reference`
- `pbm_motif`

For DBP48:

| Concept | Sequence |
|---|---|
| Original design target | `CGCCCAAAGCCGCG` |
| Experimental assay/PBM reference | `CGACACCTGACGCG` |
| PBM motif | `CTGACG` |

## Leakage graph

An edge is added between two DBPs when they share any of:

- protein sequence cluster
- original target group
- assay target group
- canonical motif group

The generated `metadata/v0_5/v0_5_split_audit.csv` is the human-readable
summary of these groups and legal split schemes. Pairwise graph edges are
preserved separately in `metadata/v0_5/v0_5_split_edges.csv`. Connected
components, rather than row counts, define the strictest known independence
unit for this seven-protein benchmark.

## Legal future splits

- Primary: leave one complete protein sequence cluster out.
- Strict sensitivity: leave one combined leakage component out.
- Illegal: random split of protein-7-mer rows.
- Illegal: split an oriented 7-mer and its reverse complement into different
  partitions.

The generated CSV manifests are the executable record of this contract.
