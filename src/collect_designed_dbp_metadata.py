from __future__ import annotations

import hashlib
import re
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import requests

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, is_valid_dna, is_valid_protein, load_yaml, normalize_sequence, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
RAW_DIR = ensure_dir(ROOT / "data" / "raw" / "gse237017")
METADATA_DIR = ensure_dir(ROOT / "metadata" / "v0_3")

SUPPLEMENTARY_TABLE_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41594-025-01669-4/"
    "MediaObjects/41594_2025_1669_MOESM3_ESM.xlsx"
)
SUPPLEMENTARY_TABLE_PATH = RAW_DIR / "nature_41594_2025_1669_MOESM3_ESM.xlsx"
PAPER_REFERENCE = "Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4; PMID 40940539"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_supplementary_table() -> Path:
    if SUPPLEMENTARY_TABLE_PATH.exists():
        return SUPPLEMENTARY_TABLE_PATH
    response = requests.get(SUPPLEMENTARY_TABLE_URL, timeout=120)
    response.raise_for_status()
    SUPPLEMENTARY_TABLE_PATH.write_bytes(response.content)
    return SUPPLEMENTARY_TABLE_PATH


def normalize_dbp_id(raw: str) -> str:
    match = re.search(r"DBP0*(\d+)", str(raw), flags=re.IGNORECASE)
    return f"DBP{match.group(1)}" if match else ""


def parse_design_table(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Supplementary Table 1", header=None)
    rows = []
    for _, row in df.iterrows():
        label = str(row.iloc[0]) if pd.notna(row.iloc[0]) else ""
        seq = str(row.iloc[1]) if pd.notna(row.iloc[1]) else ""
        if not re.search(r"DBP0*\d+", label):
            continue
        parts = [part.strip() for part in label.split(",", 1)]
        protein_id = normalize_dbp_id(parts[0])
        target_label = parts[1] if len(parts) > 1 else ""
        target_match = re.match(r"([A-Z])(?:_([bc]))?$", target_label)
        target_id = target_match.group(1) if target_match else ""
        target_context = target_match.group(2) if target_match and target_match.group(2) else ""
        protein_sequence = normalize_sequence(seq.replace("-", ""))
        rows.append(
            {
                "protein_id": protein_id,
                "design_label": label,
                "target_id": target_id,
                "target_context": target_context,
                "protein_sequence": protein_sequence,
                "sequence_length": len(protein_sequence),
                "protein_sequence_valid_20aa": is_valid_protein(protein_sequence),
            }
        )
    return pd.DataFrame(rows)


def parse_target_table(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Supplementary Table 2", header=None)
    rows = []
    for _, row in df.iterrows():
        label = str(row.iloc[0]) if pd.notna(row.iloc[0]) else ""
        sequence = str(row.iloc[1]) if pd.notna(row.iloc[1]) else ""
        if not re.match(r"^[A-Z](?:\s+\([^)]+\))?$", label.strip()):
            continue
        target_id = label.strip().split()[0]
        strands = [normalize_sequence(part) for part in sequence.split("/")]
        strand_1 = strands[0] if strands else ""
        strand_2 = strands[1] if len(strands) > 1 else ""
        rows.append(
            {
                "target_id": target_id,
                "target_label": label.strip(),
                "intended_target_dna": strand_1,
                "complementary_strand": strand_2,
                "target_duplex": f"{strand_1}/{strand_2}" if strand_2 else strand_1,
                "target_length": len(strand_1),
                "target_valid_dna": is_valid_dna(strand_1),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    path = download_supplementary_table()
    digest = sha256_file(path)
    retrieval_date = date.today().isoformat()
    selected = set(CONFIG["benchmark_v0_3"]["designed_dbp_ids"])
    designs = parse_design_table(path)
    targets = parse_target_table(path)
    selected_designs = designs[designs["protein_id"].isin(selected)].copy()
    missing = sorted(selected - set(selected_designs["protein_id"]))
    if missing:
        raise ValueError(f"Missing designed DBP sequences in supplementary table: {missing}")
    sequence_rows = []
    target_rows = []
    merged = selected_designs.merge(targets, on="target_id", how="left")
    for _, row in merged.sort_values("protein_id").iterrows():
        confidence = "high" if bool(row["protein_sequence_valid_20aa"]) else "low"
        sequence_rows.append(
            {
                "protein_id": row["protein_id"],
                "protein_sequence": row["protein_sequence"],
                "sequence_length": int(row["sequence_length"]),
                "source_type": "Nature supplementary table",
                "source_id": "Supplementary Table 1",
                "source_url": SUPPLEMENTARY_TABLE_URL,
                "paper_reference": PAPER_REFERENCE,
                "retrieval_date": retrieval_date,
                "sequence_confidence": confidence,
                "notes": f"Parsed from design label {row['design_label']}; target label {row['target_id']}_{row['target_context']}.",
            }
        )
        target_confidence = "high" if bool(row.get("target_valid_dna", False)) else "low"
        target_rows.append(
            {
                "protein_id": row["protein_id"],
                "intended_target_dna": row.get("intended_target_dna", ""),
                "target_length": int(row.get("target_length", 0)) if pd.notna(row.get("target_length", pd.NA)) else pd.NA,
                "source": "Nature supplementary table",
                "source_url": SUPPLEMENTARY_TABLE_URL,
                "paper_table_or_figure": "Supplementary Table 1 target ID joined to Supplementary Table 2 dsDNA target sequence",
                "retrieval_date": retrieval_date,
                "confidence": target_confidence,
                "notes": f"Target ID {row['target_id']} from design label {row['design_label']}; duplex={row.get('target_duplex', '')}; context suffix {row['target_context']} retained from paper label.",
                "target_id": row["target_id"],
                "target_context": row["target_context"],
                "target_duplex": row.get("target_duplex", ""),
                "complementary_strand": row.get("complementary_strand", ""),
            }
        )

    pd.DataFrame(sequence_rows).to_csv(METADATA_DIR / "designed_dbp_sequences.csv", index=False)
    pd.DataFrame(target_rows).to_csv(METADATA_DIR / "designed_dbp_targets.csv", index=False)
    manifest = pd.DataFrame(
        [
            {
                "source_type": "Nature supplementary workbook",
                "source_id": "41594_2025_1669_MOESM3_ESM.xlsx",
                "source_url": SUPPLEMENTARY_TABLE_URL,
                "local_path": str(path),
                "sha256": digest,
                "file_size_bytes": path.stat().st_size,
                "retrieval_date": retrieval_date,
                "notes": "Contains Supplementary Tables 1-3 for designed DBP sequences and target dsDNA oligos.",
            }
        ]
    )
    manifest.to_csv(METADATA_DIR / "designed_dbp_design_source_manifest.csv", index=False)
    print(pd.DataFrame(sequence_rows)[["protein_id", "sequence_length", "sequence_confidence"]].to_string(index=False))
    print(pd.DataFrame(target_rows)[["protein_id", "intended_target_dna", "target_length", "confidence"]].to_string(index=False))


if __name__ == "__main__":
    main()
