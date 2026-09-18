# v1.3.3 Structure Provenance

Retrieval date: 2026-09-18. Structure coordinates were downloaded from the
official RCSB PDB file host (`https://files.rcsb.org/download/{PDB}.cif`).
The PDB IDs and chain choices were frozen in
`V1_3_3_NATURAL_TF_CHALLENGE_CONTRACT.md` before natural-TF model performance
was inspected.

| PDB | TF candidate | Source file | Use | Primary status |
|---|---|---|---|---|
| 8OVW | Cbf1 | `data/raw/v1_3_3_external/structures/8OVW.cif` | CBF1-CCAN bound centromeric DNA | audit candidate; mapping ambiguous |
| 1AAY | Egr1/Zif268 | `data/raw/v1_3_3_external/structures/1AAY.cif` | zinc-finger DNA complex | audit candidate; mapping not evaluable |
| 2STT | ETS1 | `data/raw/v1_3_3_external/structures/2STT.cif` | ETS1/DNA NMR complex | four SaMBA sites eligible by frozen sequence mapping |
| 3WTU | ETS1 | `data/raw/v1_3_3_external/structures/3WTU.cif` | alternate ETS1/RUNX1/CBFBETA complex | retained as documented alternative, not primary |
| 1R4R | GR | `data/raw/v1_3_3_external/structures/1R4R.cif` | GR DNA-binding-domain/DNA complex | audit candidate; mapping ambiguous |
| 1HLO | Max | `data/raw/v1_3_3_external/structures/1HLO.cif` | intact human Max/DNA complex | audit candidate; mapping not evaluable |
| 1TGH | TBP | `data/raw/v1_3_3_external/structures/1TGH.cif` | human TBP/TATA DNA complex | audit candidate; mapping ambiguous |
| 1TUP | p53 | existing `data/raw/rcsb/mmcif/1TUP.cif` | p53/DNA complex | audit candidate; mapping ambiguous |

The persistent size/mtime-aware manifest records one SHA256 for each newly
downloaded file. Unchanged files are not rehashed by the v1.3.3 scripts.

Authoritative RCSB records used for candidate metadata include [1HLO](https://www.rcsb.org/structure/1HLO), [1R4R](https://www.rcsb.org/structure/1R4R), [1TGH](https://www.rcsb.org/structure/1TGH), [1AAY](https://www.rcsb.org/structure/1AAY), [2STT](https://www.rcsb.org/structure/2STT), and [8OVW](https://www.rcsb.org/structure/8OVW).
