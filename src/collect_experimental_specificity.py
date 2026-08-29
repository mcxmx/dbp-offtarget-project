from __future__ import annotations

import math
import random
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import requests

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
RAW_JASPAR_DIR = ensure_dir(ROOT / "data" / "raw" / "jaspar")
RAW_UNIPROT_DIR = ensure_dir(ROOT / "data" / "raw" / "uniprot")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")
METADATA_DIR = ensure_dir(ROOT / "metadata")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")

BASES = "ACGT"


def fetch_text(url: str, cache_path: Path) -> str:
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8")
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    cache_path.write_text(response.text, encoding="utf-8")
    return response.text


def fetch_json(url: str, cache_path: Path) -> dict[str, Any]:
    text = fetch_text(url, cache_path)
    return requests.models.complexjson.loads(text)


def fetch_jaspar_matrix(matrix_id: str) -> dict[str, Any]:
    url = f"https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/"
    return fetch_json(url, RAW_JASPAR_DIR / f"{matrix_id}.json")


def fetch_uniprot_sequence(uniprot_id: str) -> str:
    if not uniprot_id:
        return ""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    text = fetch_text(url, RAW_UNIPROT_DIR / f"{uniprot_id}.fasta")
    lines = [line.strip() for line in text.splitlines() if line and not line.startswith(">")]
    return normalize_sequence("".join(lines))


def pwm_probabilities(pfm: dict[str, list[float]], pseudocount: float) -> list[dict[str, float]]:
    length = len(pfm["A"])
    probabilities = []
    for pos in range(length):
        total = sum(float(pfm[base][pos]) for base in BASES) + pseudocount * len(BASES)
        probabilities.append({base: (float(pfm[base][pos]) + pseudocount) / total for base in BASES})
    return probabilities


def consensus_sequence(probabilities: list[dict[str, float]]) -> str:
    return "".join(max(BASES, key=lambda base: probabilities[pos][base]) for pos in range(len(probabilities)))


def pwm_log2_odds(sequence: str, probabilities: list[dict[str, float]]) -> float:
    score = 0.0
    for pos, base in enumerate(sequence):
        score += math.log2(probabilities[pos][base] / 0.25)
    return score


def sampled_sequences(consensus: str, rng: random.Random, n_random: int) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = [("consensus", consensus)]
    for pos, original in enumerate(consensus):
        for alt in BASES:
            if alt == original:
                continue
            sequence = consensus[:pos] + alt + consensus[pos + 1 :]
            rows.append(("single_substitution_from_consensus", sequence))
    length = len(consensus)
    for _ in range(n_random):
        rows.append(("random_sampled_kmer", "".join(rng.choice(BASES) for _ in range(length))))
    return rows


def source_note(matrix: dict[str, Any]) -> str:
    comments = matrix.get("comment") or []
    comment_text = "; ".join(str(item) for item in comments) if isinstance(comments, list) else str(comments)
    return (
        "JASPAR PFM-derived PWM log2-odds score; not raw PBM/HT-SELEX enrichment. "
        f"JASPAR type={matrix.get('type', '')}; comment={comment_text}"
    )


def main() -> None:
    cfg = CONFIG["experimental_specificity"]
    rng = random.Random(int(CONFIG["seed"]))
    rows = []
    source_rows = []
    for matrix_id in cfg["matrix_ids"]:
        matrix = fetch_jaspar_matrix(matrix_id)
        uniprot_ids = matrix.get("uniprot_ids") or []
        primary_uniprot = uniprot_ids[0] if uniprot_ids else ""
        protein_sequence = fetch_uniprot_sequence(primary_uniprot) if primary_uniprot else ""
        probabilities = pwm_probabilities(matrix["pfm"], float(cfg["pwm_pseudocount"]))
        consensus = consensus_sequence(probabilities)
        source_url = f"https://jaspar.elixir.no/matrix/{matrix_id}/"
        api_url = f"https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/"
        pubmed_ids = matrix.get("pubmed_ids") or []
        species = matrix.get("species") or []
        species_name = ";".join(str(item.get("name", "")) for item in species if isinstance(item, dict))
        for sequence_class, dna_sequence in sampled_sequences(
            consensus,
            rng,
            int(cfg["n_random_kmers_per_matrix"]),
        ):
            rows.append(
                {
                    "protein_id": primary_uniprot,
                    "protein_name": matrix.get("name", ""),
                    "protein_sequence": protein_sequence,
                    "dna_sequence": dna_sequence,
                    "experimental_score": pwm_log2_odds(dna_sequence, probabilities),
                    "score_type": "jaspar_pfm_pwm_log2_odds_derived",
                    "experiment_type": "motif_matrix_pfm",
                    "source_database": cfg["source_database"],
                    "source_id": matrix_id,
                    "paper_doi": "",
                    "paper_pmid": ";".join(str(pmid) for pmid in pubmed_ids),
                    "source_url": source_url,
                    "source_api_url": api_url,
                    "species": species_name,
                    "sequence_class": sequence_class,
                    "notes": source_note(matrix),
                }
            )
        source_rows.append(
            {
                "source_database": cfg["source_database"],
                "source_id": matrix_id,
                "protein_name": matrix.get("name", ""),
                "uniprot_id": primary_uniprot,
                "pubmed_ids": ";".join(str(pmid) for pmid in pubmed_ids),
                "source_url": source_url,
                "source_api_url": api_url,
                "retrieval_date": CONFIG["retrieval_date"],
                "raw_matrix_path": str(RAW_JASPAR_DIR / f"{matrix_id}.json"),
                "raw_uniprot_path": str(RAW_UNIPROT_DIR / f"{primary_uniprot}.fasta") if primary_uniprot else "",
                "notes": "Small Layer C pilot from JASPAR PFM. Scores are PWM-derived and not cross-assay normalized.",
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(PROCESSED_DIR / "experimental_specificity_small.csv", index=False)
    sources = pd.DataFrame(source_rows)
    sources.to_csv(METADATA_DIR / "experimental_specificity_sources.csv", index=False)
    summary = (
        out.groupby(["source_database", "source_id", "protein_name", "score_type"])
        .agg(
            n_sequences=("dna_sequence", "size"),
            protein_sequence_available=("protein_sequence", lambda x: bool(str(x.iloc[0]))),
            min_score=("experimental_score", "min"),
            max_score=("experimental_score", "max"),
        )
        .reset_index()
    )
    summary.to_csv(RESULTS_TABLES / "experimental_specificity_small_summary.csv", index=False)
    print(f"experimental specificity proteins: {out['protein_id'].nunique()}")
    print(f"experimental specificity rows: {len(out)}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
