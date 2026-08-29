from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import requests

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import canonical_dna, ensure_dir, gc_content, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
RETRIEVAL_DATE = CONFIG.get("retrieval_date", date.today().isoformat())

RAW_ENTRY_DIR = ensure_dir(ROOT / "data" / "raw" / "rcsb" / "entries")
RAW_ENTITY_DIR = ensure_dir(ROOT / "data" / "raw" / "rcsb" / "polymer_entities")
INTERIM_DIR = ensure_dir(ROOT / "data" / "interim")
META_DIR = ensure_dir(ROOT / "metadata")

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "dbp-offtarget-prototype/0.1"})


def fetch_json(url: str, cache_path: Path) -> dict:
    if cache_path.exists():
        with open(cache_path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    response = SESSION.get(url, timeout=60)
    response.raise_for_status()
    data = response.json()
    with open(cache_path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
    return data


def safe_join(values) -> str:
    if not values:
        return ""
    if isinstance(values, str):
        return values
    return ";".join(str(v) for v in values)


def extract_primary_citation(entry: dict) -> dict:
    citation = entry.get("rcsb_primary_citation", {}) or {}
    return {
        "source_paper_title": citation.get("title", ""),
        "source_paper_doi": citation.get("pdbx_database_id_DOI", "") or citation.get("doi", ""),
        "source_paper_pmid": citation.get("pdbx_database_id_PubMed", ""),
        "source_paper_year": citation.get("year", ""),
    }


def classify_entry(title: str) -> str:
    title = (title or "").lower()
    return "designed" if "designed" in title else "experimental"


def collect_entry(pdb_id: str) -> dict | None:
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    entry_cache = RAW_ENTRY_DIR / f"{pdb_id}.json"
    entry = fetch_json(entry_url, entry_cache)
    title = entry.get("struct", {}).get("title", "")
    poly_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
    proteins: list[dict] = []
    dnas: list[dict] = []
    for entity_id in map(str, poly_ids):
        entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
        entity_cache = RAW_ENTITY_DIR / f"{pdb_id}_{entity_id}.json"
        entity = fetch_json(entity_url, entity_cache)
        entity_poly = entity.get("entity_poly", {}) or {}
        entity_type = (entity_poly.get("type", "") or "").lower()
        seq = normalize_sequence(entity_poly.get("pdbx_seq_one_letter_code_can", ""))
        description = entity.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        chain_ids = safe_join(entity.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", []))
        record = {
            "entity_id": entity_id,
            "chain_ids": chain_ids,
            "sequence": seq,
            "sequence_length": len(seq),
            "description": description,
            "entity_type": entity_type,
            "raw_entity_json_url": entity_url,
        }
        if "polypeptide" in entity_type:
            proteins.append(record)
        elif "deoxyribonucleotide" in entity_type or "dna" in entity_type:
            dnas.append(record)
    if not proteins or not dnas:
        return None

    protein = sorted(proteins, key=lambda r: (-r["sequence_length"], r["entity_id"]))[0]
    dna = sorted(dnas, key=lambda r: (-r["sequence_length"], r["entity_id"]))[0]
    citation = extract_primary_citation(entry)

    pair_id = f"{pdb_id}_{protein['entity_id']}_{dna['entity_id']}"
    notes = (
        f"selected longest protein entity from {len(proteins)} protein entity(ies) and "
        f"longest DNA entity from {len(dnas)} DNA entity(ies)"
    )
    if "designed" in title.lower():
        notes += "; designed complex"

    row = {
        "pair_id": pair_id,
        "pdb_id": pdb_id,
        "protein_entity_id": protein["entity_id"],
        "protein_chain_ids": protein["chain_ids"],
        "protein_name": protein["description"],
        "protein_sequence": protein["sequence"],
        "protein_length": protein["sequence_length"],
        "dna_entity_id": dna["entity_id"],
        "dna_chain_ids": dna["chain_ids"],
        "target_dna": canonical_dna(dna["sequence"]),
        "dna_length": dna["sequence_length"],
        "gc_content": gc_content(dna["sequence"]),
        "source_type": "PDB",
        "source_name": "RCSB PDB",
        "source_id": pdb_id,
        "source_url": f"https://www.rcsb.org/structure/{pdb_id}",
        "source_api_url": entry_url,
        "source_title": title,
        "source_paper_title": citation["source_paper_title"],
        "source_paper_doi": citation["source_paper_doi"],
        "source_paper_pmid": citation["source_paper_pmid"],
        "source_paper_year": citation["source_paper_year"],
        "retrieval_date": RETRIEVAL_DATE,
        "experimental_or_designed": classify_entry(title),
        "has_specificity_ground_truth": True,
        "notes": notes,
    }
    return row


def main() -> None:
    rows = []
    manifest = []
    for pdb_id in CONFIG["selected_pdb_ids"]:
        try:
            row = collect_entry(pdb_id)
            if row is None:
                manifest.append({"pdb_id": pdb_id, "status": "skipped_no_dna_or_protein"})
                continue
            rows.append(row)
            manifest.append(
                {
                    "pdb_id": pdb_id,
                    "status": "collected",
                    "pair_id": row["pair_id"],
                    "source_url": row["source_url"],
                    "source_api_url": row["source_api_url"],
                    "source_title": row["source_title"],
                    "retrieval_date": RETRIEVAL_DATE,
                }
            )
        except Exception as exc:
            manifest.append({"pdb_id": pdb_id, "status": "error", "error": repr(exc)})

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["dna_length", "pair_id"], ascending=[False, True]).reset_index(drop=True)
    raw_path = INTERIM_DIR / "dbp_target_pairs_raw.csv"
    df.to_csv(raw_path, index=False)

    manifest_df = pd.DataFrame(manifest)
    manifest_df.to_csv(META_DIR / "source_manifest.csv", index=False)
    pd.DataFrame(
        [{"pdb_id": pdb_id, "selection_basis": "RCSB title search plus protein/DNA entity filter"} for pdb_id in CONFIG["selected_pdb_ids"]]
    ).to_csv(META_DIR / "selected_pdb_ids.csv", index=False)

    print(f"Collected {len(df)} protein-DNA pairs")
    print(f"Raw table written to {raw_path}")
    print(f"Source manifest written to {META_DIR / 'source_manifest.csv'}")


if __name__ == "__main__":
    main()

