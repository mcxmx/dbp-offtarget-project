# V0.6 Representation Audit

## Scope and freeze

v0.6 is a representation-repair and conditionality-audit stage. It does not
modify or overwrite `results/v0_5/`, `docs/v0_5*`, or frozen v0.5 metadata.
All seven GSE237017 designed DBPs are already `development_exposed` according
to the current exposure manifest. Therefore every real-data replay result in
this report is developmental/debugging evidence, not independent confirmatory
validation.

## Current v0.5 representation problems

v0.5 encoded each candidate by averaging the full one-hot sequence with the
full one-hot reverse complement before learning. This is RC invariant, but it
is not injective over RC classes: positional base information is averaged away.
The target encoder additionally canonicalized each sliding 7-mer, deduplicated
with a set, sorted the survivors, and mean-pooled them. That loses window order,
register, and repeated-window multiplicity.

## Exhaustive RC audit

| quantity | value |
| --- | ---: |
| oriented 7-mers | 16384 |
| canonical RC classes | 8192 |
| v0.5 representation vectors | 2000 |
| v0.5 cross-class collision groups | 1872 |
| v0.5 collision excess over RC pairs | 12384 |
| v0.6 representation vectors | 8192 |
| v0.6 cross-class collision groups | 0 |
| v0.6 RC-class injective | True |

The v0.6 candidate representation is the canonical RC representative followed
by ordinary one-hot encoding. A sequence and its reverse complement share one
vector, while different RC classes have different vectors. Independent score
invariance was tested over all 16,384 oriented 7-mers for M0, M1, M2, and M3.

## Target repair

v0.6 canonicalizes the complete target once, then retains every ordered window,
including duplicates, and appends a normalized register in `[-1, 1]`. Candidate
and target windows are compared explicitly. The default pooling is arithmetic
mean; max, log-sum-exp, softmax, and position-aware pooling are implemented as
small interpretable alternatives. No Transformer or attention stack was added.

## Synthetic conditionality gate

The synthetic benchmark contains P1/P2 and A-rich/C-rich target/candidate
preferences. All four P,T cells learned the required ranking reversal:
`P1: DNA_A > DNA_C`, `P2: DNA_C > DNA_A`, with target swaps also changing the
ranking. The synthetic gate is `GO`.

## GO / NO-GO

Representation tests: **GO**. Synthetic conditionality: **GO**. Real-data
conditionality after the frozen replay: **NO-GO** because the protein-shuffle
correlation median is `1.000000` and the target-shuffle correlation median
is `0.996871`, both compared with the pre-registered collapse threshold
`0.995`. The repair therefore fixes an information-loss bug but does not by
itself establish that the learned model uses protein/target conditionally.
