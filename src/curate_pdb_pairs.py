from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from Bio.PDB import MMCIFParser, NeighborSearch

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, gc_content, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
RAW_RCSB_DIR = ensure_dir(ROOT / "data" / "raw" / "rcsb")
RAW_MMCIF_DIR = ensure_dir(RAW_RCSB_DIR / "mmcif")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
METADATA_DIR = ensure_dir(ROOT / "metadata")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")

CONTACT_CUTOFF_ANGSTROM = 4.5


CURATION_RULES: dict[str, dict[str, Any]] = {
    "1B3T": {
        "binding_mechanism": "EBNA1 origin-binding domain bound to viral origin DNA",
        "sequence_specificity_class": "sequence_specific",
        "recommended_use": "core_benchmark",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "high",
        "curation_reason": "PDB title and primary citation identify EBNA1 bound to a viral origin DNA site.",
    },
    "1BNZ": {
        "binding_mechanism": "Sso7d archaeal chromosomal DNA-binding protein",
        "sequence_specificity_class": "non_specific",
        "recommended_use": "negative_control",
        "dna_role": "non_specific",
        "has_sequence_specificity_evidence": False,
        "curation_confidence": "high",
        "curation_reason": "PDB title/citation identify Sso7d chromosomal protein bound to DNA; not a sequence-specific TF benchmark case.",
    },
    "1D66": {
        "binding_mechanism": "GAL4 DNA-recognition domain bound to UAS-like DNA",
        "sequence_specificity_class": "sequence_specific",
        "recommended_use": "core_benchmark",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "high",
        "curation_reason": "PDB title and citation explicitly describe DNA recognition by GAL4.",
    },
    "1O3Q": {
        "binding_mechanism": "E. coli CAP/CRP-DNA recognition with indirect readout and DNA kinking",
        "sequence_specificity_class": "sequence_specific",
        "recommended_use": "core_benchmark",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "high",
        "curation_reason": "Primary citation describes CAP-DNA binding specificity and indirect readout.",
    },
    "1O3R": {
        "binding_mechanism": "E. coli CAP/CRP-DNA recognition with indirect readout and DNA kinking",
        "sequence_specificity_class": "sequence_specific",
        "recommended_use": "core_benchmark",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "high",
        "curation_reason": "Primary citation describes CAP-DNA binding specificity and indirect readout.",
    },
    "1O3S": {
        "binding_mechanism": "E. coli CAP/CRP-DNA recognition with indirect readout and DNA kinking",
        "sequence_specificity_class": "sequence_specific",
        "recommended_use": "core_benchmark",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "high",
        "curation_reason": "Primary citation describes altered CAP-DNA sequence recognition through altered DNA kinking.",
    },
    "1O3T": {
        "binding_mechanism": "E. coli CAP/CRP-DNA recognition with indirect readout and DNA kinking",
        "sequence_specificity_class": "sequence_specific",
        "recommended_use": "core_benchmark",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "high",
        "curation_reason": "Primary citation describes CAP-DNA binding specificity and indirect readout.",
    },
    "1TUP": {
        "binding_mechanism": "p53 core DNA-binding domain bound to response-element DNA",
        "sequence_specificity_class": "sequence_specific",
        "recommended_use": "core_benchmark",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "high",
        "curation_reason": "PDB title and primary citation identify a p53 tumor suppressor-DNA complex.",
    },
    "1WVL": {
        "binding_mechanism": "engineered multimeric DNA-binding protein based on Sac7d and GCN4 templates",
        "sequence_specificity_class": "designed_sequence_specific",
        "recommended_use": "auxiliary_case",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "medium",
        "curation_reason": "Primary citation describes design/characterization of a multimeric DNA-binding protein; kept auxiliary because it is engineered and not a standard natural TF case.",
    },
    "4E54": {
        "curated_protein_entity_id": "2",
        "binding_mechanism": "UV-DDB lesion-recognition complex; DDB2 directly contacts damaged DNA",
        "sequence_specificity_class": "lesion_specific",
        "recommended_use": "exclude_from_specificity_benchmark",
        "dna_role": "damaged_substrate",
        "has_sequence_specificity_evidence": False,
        "curation_confidence": "high",
        "curation_reason": "RCSB entities separate DDB1 and DDB2; chain-contact scan shows DDB2, not longest DDB1, directly contacts DNA.",
    },
    "4E5Z": {
        "curated_protein_entity_id": "2",
        "binding_mechanism": "UV-DDB lesion-recognition complex; DDB2 directly contacts damaged DNA",
        "sequence_specificity_class": "lesion_specific",
        "recommended_use": "exclude_from_specificity_benchmark",
        "dna_role": "damaged_substrate",
        "has_sequence_specificity_evidence": False,
        "curation_confidence": "high",
        "curation_reason": "RCSB entities separate DDB1 and DDB2; chain-contact scan shows DDB2, not longest DDB1, directly contacts DNA.",
    },
    "4JBM": {
        "binding_mechanism": "AIM2 HIN-domain dsDNA sensing/signaling complex",
        "sequence_specificity_class": "non_specific",
        "recommended_use": "negative_control",
        "dna_role": "non_specific",
        "has_sequence_specificity_evidence": False,
        "curation_confidence": "medium",
        "curation_reason": "Citation describes AIM2-mediated signaling termination; this is a dsDNA sensor case rather than a sequence-specific TF benchmark case.",
    },
    "4KPY": {
        "binding_mechanism": "Thermus thermophilus Argonaute guide-strand-mediated target DNA cleavage",
        "sequence_specificity_class": "guide_dependent",
        "recommended_use": "exclude_from_specificity_benchmark",
        "dna_role": "target",
        "has_sequence_specificity_evidence": False,
        "curation_confidence": "high",
        "curation_reason": "Primary citation identifies guide-strand-mediated DNA target cleavage; protein specificity is guide-dependent.",
    },
    "8TAC": {
        "binding_mechanism": "computationally designed sequence-specific DNA-binding protein",
        "sequence_specificity_class": "designed_sequence_specific",
        "recommended_use": "core_benchmark",
        "dna_role": "cognate_site",
        "has_sequence_specificity_evidence": True,
        "curation_confidence": "medium",
        "curation_reason": "PDB title/citation identify a designed sequence-specific DNA-binding protein; no quantitative specificity ground truth is attached in this repository.",
    },
    "8XIH": {
        "binding_mechanism": "Piwi/Argonaute-family protein-DNA complex with guide/target DNA components",
        "sequence_specificity_class": "guide_dependent",
        "recommended_use": "exclude_from_specificity_benchmark",
        "dna_role": "target",
        "has_sequence_specificity_evidence": False,
        "curation_confidence": "medium",
        "curation_reason": "RCSB entity annotation identifies a Piwi domain-containing protein with multiple DNA strands; treated conservatively as guide-dependent.",
    },
    "9C0F": {
        "binding_mechanism": "piggyBat transposase bound to transposon-end DNA substrate",
        "sequence_specificity_class": "structure_specific",
        "recommended_use": "method_demo_only",
        "dna_role": "substrate",
        "has_sequence_specificity_evidence": False,
        "curation_confidence": "high",
        "curation_reason": "PDB title/citation identify a transposase-substrate complex; not a standard DBP sequence-specific recognition benchmark case.",
    },
}


def read_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def clean_id(value: Any) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.endswith(".0") and text[:-2].isdigit():
        return text[:-2]
    return text


def entry_metadata(pdb_id: str) -> dict[str, str]:
    entry = read_json(RAW_RCSB_DIR / "entries" / f"{pdb_id}.json")
    citation = (entry.get("citation") or [{}])[0]
    return {
        "paper_title": citation.get("title") or "",
        "doi": citation.get("pdbx_database_id_doi") or "",
        "pmid": clean_id(citation.get("pdbx_database_id_PubMed")),
        "source_title": entry.get("struct", {}).get("title") or "",
    }


def entity_metadata(pdb_id: str, entity_id: str) -> dict[str, Any]:
    entity = read_json(RAW_RCSB_DIR / "polymer_entities" / f"{pdb_id}_{entity_id}.json")
    container = entity.get("rcsb_polymer_entity_container_identifiers", {})
    refs = container.get("reference_sequence_identifiers") or []
    uniprot_ids = [
        ref.get("database_accession", "")
        for ref in refs
        if ref.get("database_name") == "UniProt" and ref.get("database_accession")
    ]
    sequence = normalize_sequence(entity.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can"))
    return {
        "entity_id": str(entity_id),
        "entity_type": entity.get("entity_poly", {}).get("rcsb_entity_polymer_type") or "",
        "name": entity.get("rcsb_polymer_entity", {}).get("pdbx_description") or "",
        "chain_ids": ";".join(container.get("auth_asym_ids") or container.get("asym_ids") or []),
        "sequence": sequence,
        "length": len(sequence),
        "uniprot_id": ";".join(sorted(set(uniprot_ids))),
    }


def download_mmcif(pdb_id: str) -> Path:
    path = RAW_MMCIF_DIR / f"{pdb_id}.cif"
    if path.exists():
        return path
    response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=120)
    response.raise_for_status()
    with open(path, "wb") as handle:
        handle.write(response.content)
    return path


def polymer_entities(pdb_id: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    entry = read_json(RAW_RCSB_DIR / "entries" / f"{pdb_id}.json")
    proteins: list[dict[str, Any]] = []
    dnas: list[dict[str, Any]] = []
    for entity_id in entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []):
        meta = entity_metadata(pdb_id, str(entity_id))
        if meta["entity_type"] == "Protein":
            proteins.append(meta)
        elif meta["entity_type"] == "DNA":
            dnas.append(meta)
    return proteins, dnas


def contact_rows_for_entry(pdb_id: str) -> list[dict[str, Any]]:
    proteins, dnas = polymer_entities(pdb_id)
    if not proteins or not dnas:
        return []
    cif_path = download_mmcif(pdb_id)
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(cif_path))[0]
    chains = {chain.id: chain for chain in structure.get_chains()}
    rows: list[dict[str, Any]] = []
    for protein in proteins:
        for protein_chain in protein["chain_ids"].split(";"):
            if not protein_chain or protein_chain not in chains:
                continue
            protein_atoms = [
                atom
                for atom in chains[protein_chain].get_atoms()
                if (getattr(atom, "element", "") or "").upper() != "H"
            ]
            search = NeighborSearch(protein_atoms)
            for dna in dnas:
                for dna_chain in dna["chain_ids"].split(";"):
                    if not dna_chain or dna_chain not in chains:
                        continue
                    dna_atoms = [
                        atom
                        for atom in chains[dna_chain].get_atoms()
                        if (getattr(atom, "element", "") or "").upper() != "H"
                    ]
                    atom_contacts = 0
                    protein_residues = set()
                    dna_residues = set()
                    for dna_atom in dna_atoms:
                        hits = search.search(dna_atom.coord, CONTACT_CUTOFF_ANGSTROM, level="A")
                        if hits:
                            atom_contacts += len(hits)
                            dna_residues.add(dna_atom.get_parent().id)
                            for hit in hits:
                                protein_residues.add(hit.get_parent().id)
                    rows.append(
                        {
                            "pdb_id": pdb_id,
                            "protein_entity_id": protein["entity_id"],
                            "protein_chain": protein_chain,
                            "protein_name": protein["name"],
                            "dna_entity_id": dna["entity_id"],
                            "dna_chain": dna_chain,
                            "dna_name": dna["name"],
                            "atom_contacts_4p5a": atom_contacts,
                            "protein_residue_contacts_4p5a": len(protein_residues),
                            "dna_residue_contacts_4p5a": len(dna_residues),
                        }
                    )
    return rows


def best_contact(
    contact_df: pd.DataFrame,
    pdb_id: str,
    protein_entity_id: str,
    dna_entity_id: str,
) -> dict[str, Any]:
    subset = contact_df[
        (contact_df["pdb_id"] == pdb_id)
        & (contact_df["protein_entity_id"].astype(str) == str(protein_entity_id))
        & (contact_df["dna_entity_id"].astype(str) == str(dna_entity_id))
    ].copy()
    if subset.empty:
        return {"protein_chain": "", "dna_chain": "", "contacts": 0}
    subset = subset.sort_values(
        ["atom_contacts_4p5a", "protein_residue_contacts_4p5a", "dna_residue_contacts_4p5a"],
        ascending=False,
    )
    row = subset.iloc[0]
    return {
        "protein_chain": row["protein_chain"],
        "dna_chain": row["dna_chain"],
        "contacts": int(row["atom_contacts_4p5a"]),
    }


def build_curation() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")
    contact_rows: list[dict[str, Any]] = []
    for pdb_id in sorted(pairs["pdb_id"].unique()):
        contact_rows.extend(contact_rows_for_entry(str(pdb_id)))
    contact_df = pd.DataFrame(contact_rows)
    contact_df.to_csv(RESULTS_TABLES / "pdb_chain_contact_summary.csv", index=False)

    curation_rows: list[dict[str, Any]] = []
    curated_rows: list[dict[str, Any]] = []
    for _, original in pairs.iterrows():
        pdb_id = str(original["pdb_id"])
        if pdb_id not in CURATION_RULES:
            raise ValueError(f"Missing curation rule for {pdb_id}")
        rule = CURATION_RULES[pdb_id]
        original_protein_entity_id = clean_id(original["protein_entity_id"])
        original_dna_entity_id = clean_id(original["dna_entity_id"])
        curated_protein_entity_id = clean_id(rule.get("curated_protein_entity_id", original_protein_entity_id))
        curated_dna_entity_id = clean_id(rule.get("curated_dna_entity_id", original_dna_entity_id))
        protein_meta = entity_metadata(pdb_id, curated_protein_entity_id)
        dna_meta = entity_metadata(pdb_id, curated_dna_entity_id)
        paper = entry_metadata(pdb_id)
        contact = best_contact(contact_df, pdb_id, curated_protein_entity_id, curated_dna_entity_id)
        has_direct = bool(contact["contacts"] > 0)

        curation = {
            "pair_id": original["pair_id"],
            "pdb_id": pdb_id,
            "protein_name_original": original["protein_name"],
            "protein_chain_original": original["protein_chain_ids"],
            "dna_chain_original": original["dna_chain_ids"],
            "curated_protein_name": protein_meta["name"],
            "curated_protein_chain": contact["protein_chain"] or protein_meta["chain_ids"],
            "curated_dna_chain": contact["dna_chain"] or dna_meta["chain_ids"],
            "curated_protein_entity_id": curated_protein_entity_id,
            "curated_dna_entity_id": curated_dna_entity_id,
            "curated_uniprot_id": protein_meta["uniprot_id"],
            "dna_role": rule["dna_role"],
            "binding_mechanism": rule["binding_mechanism"],
            "sequence_specificity_class": rule["sequence_specificity_class"],
            "recommended_use": rule["recommended_use"],
            "has_structural_cognate": True,
            "has_direct_dna_binding_evidence": has_direct,
            "has_sequence_specificity_evidence": bool(rule["has_sequence_specificity_evidence"]),
            "has_quantitative_specificity_ground_truth": False,
            "curation_confidence": rule["curation_confidence"],
            "curation_reason": rule["curation_reason"],
            "paper_title": original.get("source_paper_title", "") or paper["paper_title"],
            "doi": original.get("source_paper_doi", "") or paper["doi"],
            "pmid": clean_id(original.get("source_paper_pmid", "")) or paper["pmid"],
            "source_url": original.get("source_url", f"https://www.rcsb.org/structure/{pdb_id}"),
            "notes": (
                f"v0.1 selected protein entity {original_protein_entity_id} and DNA entity {original_dna_entity_id}; "
                f"v0.2 curated protein entity {curated_protein_entity_id}, DNA entity {curated_dna_entity_id}; "
                f"best 4.5A contact count={contact['contacts']}. PDB structure is not quantitative specificity ground truth."
            ),
        }
        curation_rows.append(curation)

        curated = original.to_dict()
        curated.update(
            {
                "protein_entity_id": curated_protein_entity_id,
                "protein_chain_ids": curation["curated_protein_chain"],
                "protein_name": protein_meta["name"],
                "protein_sequence": protein_meta["sequence"],
                "protein_length": protein_meta["length"],
                "uniprot_id": protein_meta["uniprot_id"],
                "dna_entity_id": curated_dna_entity_id,
                "dna_chain_ids": curation["curated_dna_chain"],
                "target_dna": dna_meta["sequence"],
                "dna_length": dna_meta["length"],
                "gc_content": gc_content(dna_meta["sequence"]),
                "has_specificity_ground_truth": False,
                "has_structural_cognate": True,
                "has_direct_dna_binding_evidence": has_direct,
                "has_sequence_specificity_evidence": bool(rule["has_sequence_specificity_evidence"]),
                "has_quantitative_specificity_ground_truth": False,
                "dna_role": rule["dna_role"],
                "binding_mechanism": rule["binding_mechanism"],
                "sequence_specificity_class": rule["sequence_specificity_class"],
                "recommended_use": rule["recommended_use"],
                "curation_confidence": rule["curation_confidence"],
                "curation_reason": rule["curation_reason"],
                "notes": curation["notes"],
            }
        )
        curated_rows.append(curated)

    curation_df = pd.DataFrame(curation_rows).sort_values(["recommended_use", "pdb_id", "pair_id"])
    curated_all = pd.DataFrame(curated_rows).sort_values(["recommended_use", "pdb_id", "pair_id"])
    include_uses = CONFIG.get("benchmark_v0_2", {}).get("include_recommended_use", ["core_benchmark"])
    curated_core = curated_all[curated_all["recommended_use"].isin(include_uses)].copy()

    curation_df.to_csv(METADATA_DIR / "pdb_pair_curation.csv", index=False)
    curated_all.to_csv(PROCESSED_DIR / "dbp_target_pairs_all_curated.csv", index=False)
    curated_core.to_csv(PROCESSED_DIR / "dbp_target_pairs_v0_2.csv", index=False)
    return curation_df, curated_all, curated_core


def main() -> None:
    curation_df, curated_all, curated_core = build_curation()
    print(f"Wrote {METADATA_DIR / 'pdb_pair_curation.csv'} ({len(curation_df)} rows)")
    print(f"Wrote {PROCESSED_DIR / 'dbp_target_pairs_all_curated.csv'} ({len(curated_all)} rows)")
    print(f"Wrote {PROCESSED_DIR / 'dbp_target_pairs_v0_2.csv'} ({len(curated_core)} rows)")
    print(curated_core[["pair_id", "pdb_id", "protein_name", "sequence_specificity_class", "recommended_use"]].to_string(index=False))


if __name__ == "__main__":
    main()

