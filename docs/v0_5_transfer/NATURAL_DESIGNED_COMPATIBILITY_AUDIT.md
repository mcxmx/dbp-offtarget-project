# Natural-to-Designed Compatibility Audit

## Scope

This audit records whether the existing natural PBM benchmark can support a
new Phase 7A transfer experiment against the designed-DBP benchmark. It does
not create a new natural dataset and does not alter either frozen benchmark.

## Dataset comparison

| Dataset | Source | Proteins | Experimental units | DNA length | Primary score |
| --- | --- | ---: | ---: | ---: | --- |
| Designed | GSE237017 processed uPBM | 7 | 57,344 protein-RC-class units | 7-mer | processed PBM E-score |
| Natural | UniPROBE v0.4.1 processed benchmark | 57 | 1,875,072 protein-RC-class units | 8-mer | UniPROBE contiguous 8-mer E-score |

Both values are experimental PBM-derived enrichment-style scores, but the
datasets do not have an established identical mathematical processing
definition. The designed benchmark is a processed 7-mer E-score landscape;
the natural benchmark is a processed contiguous 8-mer UniPROBE profile.

## Compatibility findings

1. The natural benchmark is not directly assay-matched to GSE237017.
2. The natural and designed k-mer lengths differ. The project does not crop
   natural 8-mers to 7-mers, pad 7-mers to 8-mers, or invent a conversion.
3. The existing natural protein sequences are mostly full-length UniProt or
   reference sequences. Exact PBM construct/domain sequences are not
   established for the local provenance.
4. The designed benchmark contains seven curated DBP sequences and the
   independently sourced designed targets. Designed labels must not enter a
   natural-only training regime.
5. RC-class handling is compatible at the ranking-unit level, but this does
   not make the assays or score processing interchangeable.

## Prior SimplePC context

The prior v0.4.1 `SimpleProteinConditionalBaseline` artifact reports:

- natural held-out macro median Spearman: `0.3013536821293048`
- designed external macro median Spearman: `0.3616108129479222`

Those values came from a different natural-training protocol and are not a
matched Phase 7A bridge. They mix training distribution, protein
representation, objective, normalization, split, and evaluation protocol.
They cannot identify protein diversity as a causal explanation for transfer.

## Decision

Because H3 was explicitly stopped as **NOT SUPPORTED**, no new R-DESIGNED,
R-NATURAL, R-NATURAL-BUDGET-MATCHED, or R-NATURAL+DESIGNED bridge was run.
The status-only artifacts under `results/v0_5_transfer/` make this absence
explicit. No numerical transfer claim is made.

The next scientifically useful dataset work is assay/task alignment:
verified experimental constructs and, where possible, matched PBM
k-mer/score processing.
