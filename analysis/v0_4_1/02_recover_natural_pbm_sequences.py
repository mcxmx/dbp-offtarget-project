from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
import time

import pandas as pd
import requests

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root
from src.natural_pbm import EXPECTED_8MER_RC_CLASSES


ROOT = project_root()
METADATA = ensure_dir(ROOT / "metadata" / "v0_4_1")
LOGS = ensure_dir(ROOT / "logs" / "v0_4_1")
CACHE = ensure_dir(ROOT / "data" / "interim" / "v0_4_1" / "uniprot_cache_v2")


ALIAS_MAP = {
    "Bap": "bagpipe",
    "Lbl": "ladybird late",
    "Msh": "muscle segmentation homeobox",
    "Slou": "slouch",
    "Tin": "tinman",
    "Ubx": "ultrabithorax",
    "Ceh-22": "ceh-22",
    "CEH-22": "ceh-22",
    "Oct-1": "oct-1",
    "Zif268": "Zif268",
    "HLH-1": "MyoD protein 1 homolog",
    "Hmlalpha2": "MATalpha2",
    "Jun_Fos": "c-Jun",
}

SPECIES_HINTS = {
    "Busser et al., Development 2012": "Drosophila melanogaster",
    "Busser et al., PNAS 2012": "Drosophila melanogaster",
    "De Masi et al., NAR 2011": "Caenorhabditis elegans",
    "Gordan et al., Gen. Bio. 2011": "Saccharomyces cerevisiae",
    "Liu et al., Cell 2018": "Homo sapiens",
    "Liu et al., eLife 2018": "Homo sapiens",
    "Campbell et al., PLoS Pathog 2010": "Plasmodium falciparum",
    "Cheatle Jarvela et al., Mol Biol Evol 2014": "Trypanosoma brucei",
    "Del Bianco et al., PLoS ONE 2010": "Homo sapiens",
    "Alibes et al., NAR 2010": "Saccharomyces cerevisiae",
    "Berger et al., Nat Biotech 2006": "mixed",
}

EXCLUDE_PATTERNS = ["GST-", "MAML1", "NOTCH", "_rep", "_L13", "d23d456", "XL"]


@dataclass
class UniProtHit:
    accession: str
    protein_name: str
    gene_names: str
    organism: str
    length: int
    sequence: str


def cache_key(query: str) -> str:
    safe = "".join(ch if ch.isalnum() else "_" for ch in query)[:160]
    return safe + ".tsv"


def fetch_uniprot_tsv(query: str) -> str:
    cache_path = CACHE / cache_key(query)
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return cache_path.read_text(encoding="utf-8")
    last_error = None
    for attempt in range(2):
        try:
            response = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={
                    "query": query,
                    "fields": "accession,protein_name,gene_names,organism_name,length,sequence",
                    "format": "tsv",
                    "size": 10,
                },
                timeout=15,
            )
            response.raise_for_status()
            cache_path.write_text(response.text, encoding="utf-8")
            return response.text
        except requests.RequestException as exc:
            last_error = exc
            time.sleep(attempt + 1)
    raise RuntimeError(f"UniProt query failed after retries: {query}") from last_error


def query_uniprot(term: str, species: str | None = None) -> list[UniProtHit]:
    queries = []
    if species and species != "mixed":
        queries.append(f"{term} AND organism_name:{species} AND reviewed:true")
        queries.append(f"{term} AND organism_name:{species}")
    queries.append(f"{term} AND reviewed:true")
    queries.append(term)
    for query in queries:
        try:
            text = fetch_uniprot_tsv(query)
        except RuntimeError:
            continue
        lines = text.strip().splitlines()
        if len(lines) <= 1:
            continue
        rows = []
        for line in lines[1:]:
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            rows.append(
                UniProtHit(
                    accession=cols[0],
                    protein_name=cols[1],
                    gene_names=cols[2],
                    organism=cols[3],
                    length=int(cols[4]),
                    sequence=cols[5],
                )
            )
        if rows:
            return rows
    return []


def is_complex_or_fusion(name: str) -> bool:
    return any(pattern in name for pattern in EXCLUDE_PATTERNS)


def pick_hit(name: str, species: str | None, hits: list[UniProtHit]) -> UniProtHit | None:
    terms = [name, ALIAS_MAP.get(name, "")]
    terms = [term.lower().replace("_", "-") for term in terms if term]
    gene_terms = {term.replace("-", "").lower() for term in terms}
    species = (species or "").lower()
    scored: list[tuple[int, UniProtHit]] = []
    for hit in hits:
        protein_name = hit.protein_name.lower()
        gene_names = hit.gene_names.lower()
        organism = hit.organism.lower()
        has_name_evidence = any(term in protein_name or term in gene_names for term in terms)
        has_gene_evidence = any(term and term in gene_names.replace("-", "").lower() for term in gene_terms)
        if not (has_name_evidence or has_gene_evidence):
            continue
        score = 0
        if species and species != "mixed" and species in organism:
            score += 3
        if has_name_evidence:
            score += 4
        if has_gene_evidence:
            score += 4
        scored.append((score, hit))
    scored.sort(key=lambda item: (item[0], -item[1].length), reverse=True)
    if not scored or scored[0][0] <= 0:
        return None
    return scored[0][1]


def main() -> None:
    meta = pd.read_csv(METADATA / "natural_pbm_proteins.csv")
    qc_path = ROOT / "results" / "v0_4_1" / "tables" / "natural_pbm_qc_summary.csv"
    if qc_path.exists():
        qc = pd.read_csv(qc_path)
        complete_ids = set(
            qc.loc[
                (qc["n_rc_classes"] == EXPECTED_8MER_RC_CLASSES)
                & (qc["dna_qc_pass"])
                & (qc["rc_qc_pass"]),
                "natural_protein_id",
            ]
        )
    else:
        complete_ids = set(meta["natural_protein_id"])
    rows = []
    for i, (_, row) in enumerate(meta.iterrows(), start=1):
        name = str(row["protein_name"])
        publication = str(row.get("publication", ""))
        protein_id = str(row["natural_protein_id"])
        if protein_id not in complete_ids or is_complex_or_fusion(name) or "contig8mers" in name:
            rows.append(
                {
                    "natural_protein_id": protein_id,
                    "protein_name": name,
                    "protein_sequence": pd.NA,
                    "sequence_length": pd.NA,
                    "sequence_type": "unknown",
                    "sequence_match_to_assay": False,
                    "uniprot_id": pd.NA,
                    "species": pd.NA,
                    "protein_family": pd.NA,
                    "sequence_confidence": "low",
                    "source": "not_training_sequence_candidate",
                    "notes": "Skipped sequence recovery because the profile is incomplete, the label is a complex/fusion, or the protein label is not a clear factor name.",
                }
            )
            continue
        species = SPECIES_HINTS.get(publication)
        term = ALIAS_MAP.get(name, name)
        hit = pick_hit(name, species, query_uniprot(term, species))
        if hit is None and term != name:
            hit = pick_hit(name, None, query_uniprot(name, None))
        if hit is None:
            rows.append(
                {
                    "natural_protein_id": protein_id,
                    "protein_name": name,
                    "protein_sequence": pd.NA,
                    "sequence_length": pd.NA,
                    "sequence_type": row.get("sequence_type", "unknown"),
                    "sequence_match_to_assay": False,
                    "uniprot_id": pd.NA,
                    "species": pd.NA,
                    "protein_family": pd.NA,
                    "sequence_confidence": "low",
                    "source": "uniprot_unresolved",
                    "notes": "Could not resolve a conservative UniProt reference sequence.",
                }
            )
            continue
        rows.append(
            {
                "natural_protein_id": protein_id,
                "protein_name": name,
                "protein_sequence": hit.sequence,
                "sequence_length": hit.length,
                "sequence_type": "full_length",
                "sequence_match_to_assay": False,
                "uniprot_id": hit.accession,
                "species": hit.organism,
                "protein_family": pd.NA,
                "sequence_confidence": "high" if species and species.lower() in hit.organism.lower() else "medium",
                "source": f"UniProt:{term}",
                "notes": f"Reference sequence recovered from UniProt; construct-level assay sequence not claimed. UniProt protein name: {hit.protein_name}",
            }
        )
        if i % 10 == 0:
            print(f"processed {i}/{len(meta)} protein metadata rows", flush=True)
    seq = pd.DataFrame(rows)
    out = meta.drop(columns=[c for c in ["protein_sequence", "sequence_length", "sequence_type", "sequence_match_to_assay", "uniprot_id", "species", "protein_family", "sequence_confidence", "source", "notes"] if c in meta.columns], errors="ignore").merge(
        seq,
        on=["natural_protein_id", "protein_name"],
        how="left",
    )
    out.to_csv(METADATA / "natural_pbm_proteins.csv", index=False)
    seq.to_csv(METADATA / "natural_pbm_protein_sequences.csv", index=False)
    print(out[["protein_name", "sequence_type", "sequence_confidence", "uniprot_id"]].head(40).to_string(index=False))


if __name__ == "__main__":
    main()
