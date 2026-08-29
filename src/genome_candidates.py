from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, gc_content, load_yaml, normalize_sequence, project_root, reverse_complement, hamming_distance, edit_distance, sequence_identity


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
RAW_GENOME_DIR = ensure_dir(ROOT / "data" / "raw" / "genome")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed")


def download_genome(url: str, dest: Path) -> Path:
    if dest.exists():
        return dest
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()
    with open(dest, "wb") as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)
    return dest


def load_fasta_sequence(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    seq_parts = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


def find_seed_hits(genome: str, seed: str, max_hits: int) -> list[int]:
    hits = []
    start = 0
    while True:
        idx = genome.find(seed, start)
        if idx == -1:
            break
        hits.append(idx)
        if len(hits) >= max_hits:
            break
        start = idx + 1
    return hits


def retrieve_candidates(target: str, genome: str, k: int, max_candidates: int, max_seed_hits_per_seed: int) -> pd.DataFrame:
    target = normalize_sequence(target)
    rc_target = reverse_complement(target)
    L = len(target)
    seed_map = [(offset, target[offset : offset + k], "+") for offset in range(0, max(1, L - k + 1))]
    seed_map += [(offset, rc_target[offset : offset + k], "-") for offset in range(0, max(1, L - k + 1))]

    candidate_support = defaultdict(lambda: {"support_count": 0, "support_strands": set(), "seed_examples": []})
    for offset, seed, strand in seed_map:
        if len(seed) < k:
            continue
        for hit in find_seed_hits(genome, seed, max_seed_hits_per_seed):
            start = hit - offset
            if start < 0 or start + L > len(genome):
                continue
            info = candidate_support[start]
            info["support_count"] += 1
            info["support_strands"].add(strand)
            if len(info["seed_examples"]) < 3:
                info["seed_examples"].append(f"{strand}:{offset}:{seed}")

    rows = []
    for start, info in candidate_support.items():
        candidate = genome[start : start + L]
        rev_candidate = reverse_complement(candidate)
        h_fwd = hamming_distance(target, candidate)
        h_rev = hamming_distance(target, rev_candidate)
        if h_fwd <= h_rev:
            strand = "+"
            best_candidate = candidate
            best_hamming = h_fwd
        else:
            strand = "-"
            best_candidate = rev_candidate
            best_hamming = h_rev
        rows.append(
            {
                "chromosome": CONFIG["genome_demo"]["chromosome"],
                "start": start + 1,
                "end": start + L,
                "strand": strand,
                "genomic_sequence": candidate,
                "candidate_sequence": best_candidate,
                "hamming_distance": best_hamming,
                "edit_distance": edit_distance(target, best_candidate),
                "sequence_identity": sequence_identity(target, best_candidate),
                "gc_content": gc_content(candidate),
                "support_count": info["support_count"],
                "support_strands": ";".join(sorted(info["support_strands"])),
                "seed_examples": "|".join(info["seed_examples"]),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df.sort_values(["hamming_distance", "support_count", "start"], ascending=[True, False, True]).reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df.head(max_candidates)


def main() -> None:
    genome_path = RAW_GENOME_DIR / "chr22.fa.gz"
    download_genome(CONFIG["genome_demo"]["fasta_url"], genome_path)
    genome = load_fasta_sequence(genome_path)
    pairs = pd.read_csv(PROCESSED_DIR / "dbp_target_pairs.csv")
    row = pairs.loc[pairs["pair_id"] == CONFIG["genome_scan_pair_id"]]
    if row.empty:
        raise ValueError(f"Genome scan pair_id not found: {CONFIG['genome_scan_pair_id']}")
    row = row.iloc[0]
    target = normalize_sequence(row["target_dna"])
    candidates = retrieve_candidates(
        target,
        genome,
        int(CONFIG["genome_demo"]["seed_kmer_length"]),
        int(CONFIG["genome_demo"]["max_candidates"]),
        int(CONFIG["genome_demo"]["max_seed_hits_per_seed"]),
    )
    candidates.insert(0, "pair_id", row["pair_id"])
    candidates.insert(1, "target_dna", target)
    candidates.insert(2, "protein_name", row["protein_name"])
    candidates.insert(3, "protein_sequence", row["protein_sequence"])
    candidates["source"] = f"{CONFIG['genome_demo']['assembly']} {CONFIG['genome_demo']['chromosome']} candidate retrieval demo"
    candidates["candidate_type"] = "genome_candidate"
    candidates.to_csv(RESULTS_TABLES / "genome_candidates_demo.csv", index=False)
    print(f"Genome candidates written to {RESULTS_TABLES / 'genome_candidates_demo.csv'}")
    print(candidates.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
