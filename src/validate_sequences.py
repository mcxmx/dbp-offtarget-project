from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import (
    canonical_dna,
    ensure_dir,
    is_valid_dna,
    is_valid_protein,
    normalize_sequence,
    project_root,
    reverse_complement,
    reverse_complement_canonical,
)


ROOT = project_root()
INTERIM_DIR = ensure_dir(ROOT / "data" / "interim")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "tables")


def row_flag(reason: list[str]) -> str:
    return ";".join(reason) if reason else ""


def main() -> None:
    raw_path = INTERIM_DIR / "dbp_target_pairs_raw.csv"
    if not raw_path.exists():
        raise FileNotFoundError(f"Missing raw table: {raw_path}")

    df = pd.read_csv(raw_path)
    if df.empty:
        raise ValueError("Raw pair table is empty")

    df["protein_sequence"] = df["protein_sequence"].map(normalize_sequence)
    df["target_dna"] = df["target_dna"].map(normalize_sequence)
    df["dna_length"] = df["target_dna"].map(len)
    df["protein_length"] = df["protein_sequence"].map(len)
    df["dna_is_valid"] = df["target_dna"].map(is_valid_dna)
    df["protein_is_valid"] = df["protein_sequence"].map(is_valid_protein)
    df["dna_has_ambiguous_bases"] = ~df["dna_is_valid"]
    df["protein_has_invalid_chars"] = ~df["protein_is_valid"]
    df["dna_missing"] = df["target_dna"].eq("")
    df["protein_missing"] = df["protein_sequence"].eq("")

    df["dna_canonical"] = df["target_dna"].map(canonical_dna)
    df["dna_rc_canonical"] = df["target_dna"].map(reverse_complement_canonical)
    df["dna_reverse_complement"] = df["target_dna"].map(reverse_complement)
    df["dna_exact_duplicate"] = df["target_dna"].duplicated(keep=False)
    df["dna_rc_duplicate"] = df["dna_rc_canonical"].duplicated(keep=False)
    df["protein_exact_duplicate"] = df["protein_sequence"].duplicated(keep=False)
    df["dna_unique_group_size"] = df.groupby("dna_canonical")["pair_id"].transform("count")
    df["protein_unique_group_size"] = df.groupby("protein_sequence")["pair_id"].transform("count")

    reasons = []
    for _, row in df.iterrows():
        row_reasons = []
        if row["dna_missing"]:
            row_reasons.append("missing_dna")
        if row["protein_missing"]:
            row_reasons.append("missing_protein")
        if not row["dna_is_valid"] and not row["dna_missing"]:
            row_reasons.append("invalid_dna")
        if not row["protein_is_valid"] and not row["protein_missing"]:
            row_reasons.append("invalid_protein")
        reasons.append(row_flag(row_reasons))
    df["filter_reason"] = reasons
    df["keep_row"] = df["filter_reason"].eq("")

    validated_path = INTERIM_DIR / "dbp_target_pairs_validated.csv"
    df.to_csv(validated_path, index=False)

    summary_rows = [
        {"metric": "raw_rows", "value": int(len(df))},
        {"metric": "valid_rows", "value": int(df["keep_row"].sum())},
        {"metric": "filtered_rows", "value": int((~df["keep_row"]).sum())},
        {"metric": "invalid_dna_rows", "value": int((~df["dna_is_valid"]).sum())},
        {"metric": "invalid_protein_rows", "value": int((~df["protein_is_valid"]).sum())},
        {"metric": "missing_dna_rows", "value": int(df["dna_missing"].sum())},
        {"metric": "missing_protein_rows", "value": int(df["protein_missing"].sum())},
        {"metric": "exact_duplicate_dna_rows", "value": int(df["dna_exact_duplicate"].sum())},
        {"metric": "reverse_complement_duplicate_dna_rows", "value": int(df["dna_rc_duplicate"].sum())},
        {"metric": "exact_duplicate_protein_rows", "value": int(df["protein_exact_duplicate"].sum())},
        {"metric": "unique_dna_canonical_sequences", "value": int(df["dna_canonical"].nunique())},
        {"metric": "unique_protein_sequences", "value": int(df["protein_sequence"].nunique())},
    ]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(RESULTS_TABLES / "dataset_qc.csv", index=False)

    df[[
        "pair_id",
        "pdb_id",
        "dna_length",
        "protein_length",
        "dna_is_valid",
        "protein_is_valid",
        "dna_exact_duplicate",
        "dna_rc_duplicate",
        "protein_exact_duplicate",
        "filter_reason",
        "keep_row",
    ]].to_csv(RESULTS_TABLES / "dataset_qc_rows.csv", index=False)

    print(summary.to_string(index=False))
    print(f"Validated table written to {validated_path}")


if __name__ == "__main__":
    main()
