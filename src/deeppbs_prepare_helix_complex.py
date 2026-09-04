"""Prepare a DeepPBS input using the DSSR-defined DNA helix only.

This is a project-side input preparation step for structures containing DNA
overhangs or non-helical nucleotides. The selected DNA residues come directly
from the first DSSR helix's base-pair records; no sequence or coordinates are
invented. Protein residues are retained from the original structure.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select


def _parse_nt_id(nt_id: str) -> tuple[str, tuple[str, int, str]]:
    fields = nt_id.split(".")
    chain_id = fields[2]
    residue_id = (" ", int(fields[4]), fields[5] or " ")
    return chain_id, residue_id


class _HelixComplexSelect(Select):
    def __init__(self, protein_chains: set[str], dna_residues: set[tuple[str, tuple[str, int, str]]]):
        self.protein_chains = protein_chains
        self.dna_residues = dna_residues

    def accept_chain(self, chain):
        return chain.id in self.protein_chains or any(cid == chain.id for cid, _ in self.dna_residues)

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        residue_id = residue.id
        resname = residue.get_resname().strip()
        if chain_id in self.protein_chains and resname not in {"HOH", "WAT"}:
            return 1
        return int((chain_id, residue_id) in self.dna_residues)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("structure")
    parser.add_argument("dssr_json")
    parser.add_argument("output")
    parser.add_argument("--protein-chain", action="append", required=True)
    args = parser.parse_args()

    structure_path = Path(args.structure)
    parser_cls = MMCIFParser if structure_path.suffix.lower() == ".cif" else PDBParser
    structure = parser_cls(QUIET=True).get_structure("helix_complex", str(structure_path))
    with open(args.dssr_json, "r", encoding="utf-8") as handle:
        dssr = json.load(handle)

    helices = dssr.get("helices", [])
    if len(helices) != 1:
        raise ValueError(f"Expected exactly one DSSR helix, found {len(helices)}")
    pairs = helices[0].get("pairs", [])
    dna_residues: set[tuple[str, tuple[str, int, str]]] = set()
    for pair in pairs:
        for key in ("nt1", "nt2"):
            dna_residues.add(_parse_nt_id(pair[key]))

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), _HelixComplexSelect(set(args.protein_chain), dna_residues))

    print(f"selected_base_pairs={len(pairs)}")
    print(f"selected_dna_residues={len(dna_residues)}")
    print(f"protein_chains={','.join(args.protein_chain)}")
    print(f"output={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
