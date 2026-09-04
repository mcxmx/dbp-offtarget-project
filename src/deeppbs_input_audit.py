"""Stage-by-stage diagnostic for DeepPBS structure preprocessing.

This module is run in the upstream DeepPBS environment against a copied input
structure. It reports the first failing preprocessing stage without changing
the upstream implementation or imputing missing coordinates.
"""

from __future__ import annotations

import argparse
import inspect
import os
import json
import subprocess
import sys
import traceback
from typing import Any, Callable


def _stage(name: str, function: Callable[[], Any]) -> Any:
    print(f"\n=== {name} ===", flush=True)
    try:
        value = function()
    except Exception as exc:
        print(f"{name}: FAILED: {type(exc).__name__}: {exc}", flush=True)
        traceback.print_exc()
        raise
    print(f"{name}: OK", flush=True)
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("structure_path")
    parser.add_argument("--no-cleanp", action="store_true")
    parser.add_argument("--dump-dir")
    args = parser.parse_args()

    upstream_root = os.environ.get("DEEPPBS_ROOT", "/home/qwqaq/DeepPBS")
    sys.path.insert(0, upstream_root)

    from deeppbs import (  # type: ignore
        StructureData,
        cleanDNA,
        cleanProtein,
        countContacts,
        getAchtleyFactors,
        getAtomSASA,
        getCV,
        makeDNACG,
        makeProteinGraph,
        processDNA,
        splitEntities,
    )
    print(f"processDNA_signature={inspect.signature(processDNA)}", flush=True)
    print(f"makeDNACG_signature={inspect.signature(makeDNACG)}", flush=True)
    print(f"makeProteinGraph_signature={inspect.signature(makeProteinGraph)}", flush=True)
    print(f"countContacts_signature={inspect.signature(countContacts)}", flush=True)

    structure = _stage(
        "StructureData",
        lambda: StructureData(args.structure_path),
    )
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            dna_map = {"DA": "A", "DC": "C", "DG": "G", "DT": "T"}
            chain_sequence = "".join(
                dna_map.get(residue.get_resname().strip(), "?")
                for residue in residues
            )
            print(
                f"chain_summary model={model.id} chain={chain.id} "
                f"n_residues={len(residues)} sequence={chain_sequence} "
                f"residue_names={sorted(set(residue.get_resname().strip() for residue in residues))}",
                flush=True,
            )
    protein, dna = _stage(
        "splitEntities",
        lambda: splitEntities(structure),
    )
    print(
        f"protein_type={type(protein).__name__} "
        f"dna_type={type(dna).__name__}",
        flush=True,
    )
    for label, value in (("protein", protein), ("dna", dna)):
        public_attrs = [
            name
            for name in dir(value)
            if not name.startswith("_") and not callable(getattr(value, name, None))
        ]
        print(f"{label}_attrs={public_attrs}", flush=True)

    if args.no_cleanp:
        cleaned_protein = protein
        print("cleanProtein: SKIPPED (--no-cleanp)", flush=True)
    else:
        cleaned_protein_result = _stage("cleanProtein", lambda: cleanProtein(protein))
        if isinstance(cleaned_protein_result, tuple):
            cleaned_protein = cleaned_protein_result[0]
            print(
                "cleanProtein_return_tuple="
                f"{[type(item).__name__ for item in cleaned_protein_result]}",
                flush=True,
            )
        else:
            cleaned_protein = cleaned_protein_result
    cleaned_dna = _stage("cleanDNA", lambda: cleanDNA(dna))
    if args.dump_dir:
        os.makedirs(args.dump_dir, exist_ok=True)
        cleaned_dna_path = os.path.join(args.dump_dir, "cleaned_dna.pdb")
        cleaned_dna.save(cleaned_dna_path)
        print(f"cleaned_dna_path={cleaned_dna_path}", flush=True)
        dssr_path = os.path.join(args.dump_dir, "cleaned_dna_dssr.json")
        dssr_exe = "x3dna-dssr"
        command = [
            dssr_exe,
            f"--i={cleaned_dna_path}",
            f"--o={dssr_path}",
            "--json",
            "--more",
            "--idstr=long",
            "--non-pair",
        ]
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        print(f"manual_dssr_returncode={completed.returncode}", flush=True)
        if completed.stdout:
            print("manual_dssr_stdout=" + completed.stdout.strip(), flush=True)
        if completed.stderr:
            print("manual_dssr_stderr=" + completed.stderr.strip(), flush=True)
        if os.path.exists(dssr_path):
            with open(dssr_path, "r", encoding="utf-8") as handle:
                dssr = json.load(handle)
            missing_frame = []
            for nucleotide in dssr.get("nts", []):
                if "frame" not in nucleotide:
                    missing_frame.append(
                        {
                            "nt_id": nucleotide.get("nt_id"),
                            "nt_name": nucleotide.get("nt_name"),
                            "chain_name": nucleotide.get("chain_name"),
                            "keys": sorted(nucleotide),
                        }
                    )
            print(f"manual_dssr_nt_count={len(dssr.get('nts', []))}", flush=True)
            print(f"manual_dssr_missing_frame={missing_frame}", flush=True)
    dna_data = _stage("processDNA", lambda: processDNA(cleaned_dna, quiet=False))
    print(f"dna_data_type={type(dna_data).__name__}", flush=True)
    print(f"dna_data_len={len(dna_data)}", flush=True)
    if dna_data:
        print(f"dna_data_item_type={type(dna_data[0]).__name__}", flush=True)
        if isinstance(dna_data[0], dict):
            print(f"dna_data_item_keys={sorted(dna_data[0])}", flush=True)
            entities = dna_data[0].get("entities", [])
            helices = []
            for entity in entities:
                print(
                    f"entity_keys={sorted(entity)} "
                    f"helical_segments={len(entity.get('helical_segments', []))}",
                    flush=True,
                )
                helices.extend(entity.get("helical_segments", []))
            print(f"n_helices={len(helices)}", flush=True)
            if helices:
                print(f"helix_keys={sorted(helices[0])}", flush=True)
        else:
            helices = []
    else:
        helices = []

    dna_cg = _stage(
        "makeDNACG",
        lambda: makeDNACG(cleaned_dna, helices[0]),
    )
    helix = helices[0]
    print(
        "shape_parameter_lengths="
        + str(
            {
                key: len(value)
                for key, value in helix["shape_parameters"].items()
            }
        ),
        flush=True,
    )
    dna_feature_names = []
    feature_matrix = []
    for index in range(helix["length"]):
        row = []
        for parameter in (
            "buckle",
            "shear",
            "stretch",
            "stagger",
            "propeller",
            "opening",
        ):
            row.append(helix["shape_parameters"][parameter][index])
            dna_feature_names.append(parameter)
        feature_matrix.append(row)
    print(
        "basic_shape_matrix="
        f"shape=({len(feature_matrix)},{len(feature_matrix[0]) if feature_matrix else 0})",
        flush=True,
    )
    protein_features = ["charge", "radius"]
    sasa_feature = _stage(
        "getAtomSASA",
        lambda: getAtomSASA(cleaned_protein, classifier=None),
    )
    protein_features.append(sasa_feature)
    achtley_features = _stage(
        "getAchtleyFactors",
        lambda: getAchtleyFactors(cleaned_protein),
    )
    protein_features += achtley_features
    cv_feature = _stage(
        "getCV",
        lambda: getCV(cleaned_protein, 7.5, feature_name="cv", impute_hydrogens=True),
    )
    protein_features.append(cv_feature)
    print(f"protein_features={protein_features}", flush=True)
    dna_graph = _stage(
        "makeProteinGraph",
        lambda: makeProteinGraph(
            cleaned_protein,
            feature_names=protein_features,
        ),
    )
    contacts = _stage(
        "countContacts",
        lambda: countContacts(cleaned_protein, args.structure_path, dna_cg[0], [True] * dna_cg[0].shape[0]),
    )
    print(
        "summary="
        f"dna_cg_type={type(dna_cg).__name__},"
        f"dna_graph_type={type(dna_graph).__name__},"
        f"contacts={contacts}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
