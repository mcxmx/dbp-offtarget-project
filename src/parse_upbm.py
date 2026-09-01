from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, is_valid_dna, load_yaml, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
EXPECTED_7MERS = int(CONFIG["benchmark_v0_3"]["expected_7mer_count"])
RAW_DIR = ensure_dir(ROOT / "data" / "raw" / "gse237017")
INTERIM_DIR = ensure_dir(ROOT / "data" / "interim" / "gse237017")
METADATA_DIR = ensure_dir(ROOT / "metadata" / "v0_3")
RESULTS_TABLES = ensure_dir(ROOT / "results" / "v0_3" / "tables")


def normalize_header(column: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", str(column).strip().lower()).strip("_")
    return text


def read_7mer_file(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        df = pd.read_csv(handle, sep="\t")
    columns = list(df.columns)
    normalized = [normalize_header(col) for col in columns]
    rename = {}
    seen_7mer = 0
    for original, norm in zip(columns, normalized):
        if norm == "7_mer" or norm.startswith("7_mer_"):
            seen_7mer += 1
            if seen_7mer == 1:
                rename[original] = "dna_7mer"
            elif seen_7mer == 2:
                rename[original] = "dna_7mer_reverse_complement_column"
            else:
                rename[original] = f"dna_7mer_extra_{seen_7mer}"
        elif norm in {"e_score", "escore"}:
            rename[original] = "e_score"
        elif norm == "median":
            rename[original] = "median_intensity"
        elif norm == "z_score":
            rename[original] = "z_score"
        else:
            rename[original] = norm
    df = df.rename(columns=rename)
    for required in ["dna_7mer", "e_score", "median_intensity", "z_score"]:
        if required not in df.columns:
            df[required] = pd.NA
    df["dna_7mer"] = df["dna_7mer"].astype(str).str.upper()
    if "dna_7mer_reverse_complement_column" in df.columns:
        df["dna_7mer_reverse_complement_column"] = df["dna_7mer_reverse_complement_column"].astype(str).str.upper()
    for score_col in ["e_score", "median_intensity", "z_score"]:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df.insert(0, "source_row_index", range(1, len(df) + 1))
    return df


def expand_reverse_complement_columns(df: pd.DataFrame) -> pd.DataFrame:
    if "dna_7mer_reverse_complement_column" not in df.columns:
        out = df.copy()
        out["source_orientation"] = "primary_column"
        out["paired_7mer"] = pd.NA
        return out
    primary = df.copy()
    primary["source_orientation"] = "primary_column"
    primary["paired_7mer"] = primary["dna_7mer_reverse_complement_column"]
    rc = df.copy()
    rc["paired_7mer"] = rc["dna_7mer"]
    rc["dna_7mer"] = rc["dna_7mer_reverse_complement_column"]
    rc["source_orientation"] = "reverse_complement_column"
    return pd.concat([primary, rc], ignore_index=True, sort=False)


def parse_processed_files() -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = pd.read_csv(METADATA_DIR / "gse237017_samples.csv", dtype={"replicate": str})
    manifest = pd.read_csv(METADATA_DIR / "gse237017_file_manifest.csv")
    processed = manifest[manifest["file_type"] == "processed_7mer"].copy()
    rows = []
    qc_rows = []
    for _, file_row in processed.iterrows():
        local_path = Path(file_row["local_path"])
        if not local_path.exists():
            raise FileNotFoundError(f"Missing processed 7-mer file: {local_path}")
        sample = samples.loc[samples["gsm_id"] == file_row["gsm_id"]]
        if sample.empty:
            raise ValueError(f"No sample metadata for {file_row['gsm_id']}")
        sample_row = sample.iloc[0]
        source_df = read_7mer_file(local_path)
        source_valid = source_df["dna_7mer"].map(lambda seq: len(seq) == 7 and is_valid_dna(seq))
        if "dna_7mer_reverse_complement_column" in source_df.columns:
            source_valid = source_valid & source_df["dna_7mer_reverse_complement_column"].map(lambda seq: len(seq) == 7 and is_valid_dna(seq))
        df = expand_reverse_complement_columns(source_df)
        valid_dna = df["dna_7mer"].map(lambda seq: len(seq) == 7 and is_valid_dna(seq))
        bad_rows = df.loc[~valid_dna].copy()
        df = df.assign(
            gse_id=sample_row["gse_id"],
            gsm_id=sample_row["gsm_id"],
            protein_id=sample_row["protein_id"],
            protein_name=sample_row["protein_name"],
            protein_concentration=sample_row["protein_concentration"],
            replicate=str(sample_row["replicate"]),
            sample_title=sample_row["sample_title"],
            platform_id=sample_row["platform_id"],
            source_file=file_row["filename"],
            source_file_sha256=file_row["sha256"],
            source_url=sample_row["source_url"],
        )
        rows.append(
            df[
                [
                    "gse_id",
                    "gsm_id",
                    "protein_id",
                    "protein_name",
                    "protein_concentration",
                    "replicate",
                    "sample_title",
                    "platform_id",
                    "source_row_index",
                    "dna_7mer",
                    "dna_7mer_reverse_complement_column",
                    "paired_7mer",
                    "source_orientation",
                    "e_score",
                    "median_intensity",
                    "z_score",
                    "source_file",
                    "source_file_sha256",
                    "source_url",
                ]
            ]
        )
        n_rows = len(df)
        n_unique = int(df["dna_7mer"].nunique())
        n_source_rows = len(source_df)
        n_unique_primary = int(source_df["dna_7mer"].nunique())
        n_unique_source_union = int(
            pd.concat(
                [
                    source_df["dna_7mer"],
                    source_df["dna_7mer_reverse_complement_column"]
                    if "dna_7mer_reverse_complement_column" in source_df.columns
                    else pd.Series(dtype=str),
                ],
                ignore_index=True,
            ).nunique()
        )
        qc_rows.append(
            {
                "protein_id": sample_row["protein_id"],
                "gsm_id": sample_row["gsm_id"],
                "replicate": str(sample_row["replicate"]),
                "n_rows": n_rows,
                "n_unique_7mers": n_unique,
                "n_missing_7mers": max(EXPECTED_7MERS - n_unique, 0),
                "n_duplicate_7mers": n_rows - n_unique,
                "fraction_complete": n_unique / EXPECTED_7MERS,
                "n_source_rows": n_source_rows,
                "n_unique_primary_column_7mers": n_unique_primary,
                "n_unique_source_union_7mers": n_unique_source_union,
                "n_invalid_dna_rows": int(len(bad_rows)),
                "n_invalid_source_rows": int((~source_valid).sum()),
                "n_missing_e_score": int(df["e_score"].isna().sum()),
                "n_missing_median_intensity": int(df["median_intensity"].isna().sum()),
                "n_missing_z_score": int(df["z_score"].isna().sum()),
                "source_file": file_row["filename"],
                "notes": "GEO 7-mer table has two 7-mer columns; rows were explicitly expanded so reverse-complement companion sequences are represented as separate 7-mer measurements with the same PBM scores.",
            }
        )
    long_df = pd.concat(rows, ignore_index=True)
    qc_df = pd.DataFrame(qc_rows).sort_values(["protein_id", "replicate", "gsm_id"])
    return long_df, qc_df


def main() -> None:
    long_df, qc_df = parse_processed_files()
    long_df.to_parquet(INTERIM_DIR / "upbm_7mers_long.parquet", index=False)
    long_df.head(200).to_csv(INTERIM_DIR / "upbm_7mers_long_preview.csv", index=False)
    qc_df.to_csv(RESULTS_TABLES / "sample_coverage_qc.csv", index=False)
    print(f"parsed rows: {len(long_df)}")
    print(f"samples: {long_df['gsm_id'].nunique()}")
    print(f"proteins: {long_df['protein_id'].nunique()}")
    print(qc_df.to_string(index=False))


if __name__ == "__main__":
    main()
