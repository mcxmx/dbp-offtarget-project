from __future__ import annotations

import hashlib
import io
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.utils import DNA_ALPHABET, ensure_dir, normalize_sequence, project_root


UNIPROBE_DOWNLOADS_URL = "https://thebrain.bwh.harvard.edu/uniprobe/downloads.php"
UNIPROBE_BASE_URL = "https://thebrain.bwh.harvard.edu/uniprobe/"
EXPECTED_8MER_RC_CLASSES = (4**8 + 4**4) // 2


@dataclass(frozen=True)
class UniProbePublication:
    publication: str
    dataset_code: str
    contig8mer_url: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_uniprobe_downloads_page(html: str) -> list[UniProbePublication]:
    pattern = re.compile(
        r'<h3 class="pub_a"><a class="not-link-down" href="">(?P<pub>.*?)</a></h3>.*?'
        r"(?P<link>downloads/(?P<code>[^/]+)/(?P=code)_contig8mers\.zip)",
        flags=re.DOTALL,
    )
    publications = []
    for match in pattern.finditer(html):
        publication = re.sub(r"\s+", " ", match.group("pub")).strip()
        code = match.group("code")
        link = urljoin(UNIPROBE_BASE_URL, match.group("link"))
        publications.append(UniProbePublication(publication=publication, dataset_code=code, contig8mer_url=link))
    return publications


def download_file(url: str, out_path: Path, timeout: int = 120) -> str:
    ensure_dir(out_path.parent)
    if out_path.exists() and out_path.stat().st_size > 0:
        return "exists"
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with out_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return "downloaded"


def zip_entries(zip_path: Path) -> list[zipfile.ZipInfo]:
    with zipfile.ZipFile(zip_path) as archive:
        return [entry for entry in archive.infolist() if not entry.is_dir()]


def infer_uniprobe_ids(entry_name: str, dataset_code: str) -> dict[str, str]:
    parts = Path(entry_name).parts
    factor_name = parts[0] if parts else Path(entry_name).stem
    stem = Path(entry_name).stem
    experiment_id = re.sub(r"_8mers_11111111$|_contig8mers_combined$|_contig8mers$|_8mers$", "", stem)
    replicate_match = re.search(r"(?:^|_)rep(?:licate)?[_-]?(\d+)(?:_|$)", experiment_id, flags=re.IGNORECASE)
    replicate_id = f"rep{replicate_match.group(1)}" if replicate_match else "single_measurement"
    construct = re.sub(r"(?:_)?rep(?:licate)?[_-]?\d+(?:_|$)", "", experiment_id, flags=re.IGNORECASE).strip("_")
    natural_protein_id = f"UniPROBE:{dataset_code}:{construct}"
    return {
        "natural_protein_id": natural_protein_id,
        "protein_name": factor_name,
        "experiment_id": f"{dataset_code}:{experiment_id}",
        "replicate_id": replicate_id,
        "construct_label": construct,
    }


def read_uniprobe_8mer_file(raw_bytes: bytes) -> pd.DataFrame:
    text = raw_bytes.decode("utf-8", errors="replace").strip()
    first_line = text.splitlines()[0] if text else ""
    has_header = "E-score" in first_line or first_line.startswith("8-mer")
    df = pd.read_csv(io.StringIO(text), sep="\t", header=0 if has_header else None)
    if df.shape[1] == 3:
        df.columns = ["8-mer", "8-mer.1", "E-score"]
    elif df.shape[1] == 4:
        df.columns = ["8-mer", "8-mer.1", "E-score", "Median"]
    elif df.shape[1] >= 5:
        df.columns = ["8-mer", "8-mer.1", "E-score", "Median", "Z-score"] + [f"extra_{i}" for i in range(df.shape[1] - 5)]
    else:
        raise ValueError(f"Unexpected UniPROBE 8-mer column count: {df.shape[1]}")
    rc_col = "8-mer.1" if "8-mer.1" in df.columns else None
    out = pd.DataFrame(
        {
            "dna_sequence": df["8-mer"].astype(str).map(normalize_sequence),
            "reverse_complement_from_source": df[rc_col].astype(str).map(normalize_sequence) if rc_col else pd.NA,
            "experimental_score": pd.to_numeric(df["E-score"], errors="coerce"),
            "median_intensity": pd.to_numeric(df["Median"], errors="coerce") if "Median" in df.columns else np.nan,
            "z_score": pd.to_numeric(df["Z-score"], errors="coerce") if "Z-score" in df.columns else np.nan,
        }
    )
    out["dna_length"] = out["dna_sequence"].str.len()
    out["canonical_rc"] = out["dna_sequence"].map(canonical_rc)
    return out


def is_candidate_uniprobe_8mer_entry(name: str) -> bool:
    lower = name.lower()
    return lower.endswith(".txt") and "8mer" in lower and "top_enrichment" not in lower


def parse_uniprobe_zip(zip_path: Path, dataset_code: str, publication: str, max_profiles: int | None = None) -> pd.DataFrame:
    rows = []
    with zipfile.ZipFile(zip_path) as archive:
        names = [name for name in archive.namelist() if is_candidate_uniprobe_8mer_entry(name)]
        if max_profiles is not None:
            names = names[:max_profiles]
        for name in names:
            ids = infer_uniprobe_ids(name, dataset_code)
            table = read_uniprobe_8mer_file(archive.read(name))
            table["natural_protein_id"] = ids["natural_protein_id"]
            table["protein_name"] = ids["protein_name"]
            table["construct_label"] = ids["construct_label"]
            table["experiment_id"] = ids["experiment_id"]
            table["replicate_id"] = ids["replicate_id"]
            table["publication"] = publication
            table["database"] = "UniPROBE"
            table["assay_type"] = "universal protein-binding microarray"
            table["experimental_score_type"] = "UniPROBE contiguous 8-mer E-score"
            table["source_file"] = name
            rows.append(table)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def parse_uniprobe_archives(manifest: pd.DataFrame, root: Path | None = None) -> pd.DataFrame:
    root = root or project_root()
    rows = []
    for _, entry in manifest.iterrows():
        if entry.get("download_status") not in {"downloaded", "exists"}:
            continue
        local_path = root / str(entry["local_relative_path"])
        if not local_path.exists():
            continue
        df = parse_uniprobe_zip(local_path, str(entry["dataset_code"]), str(entry["publication"]))
        if not df.empty:
            df["archive_dataset_code"] = entry["dataset_code"]
            df["archive_publication"] = entry["publication"]
            rows.append(df)
    if not rows:
        return pd.DataFrame()
    long_df = pd.concat(rows, ignore_index=True)
    long_df["dna_valid"] = long_df["dna_sequence"].str.fullmatch("[ACGT]{8}")
    long_df["rc_source_matches"] = long_df["reverse_complement_from_source"].eq(long_df["dna_sequence"].map(reverse_complement))
    return long_df


def validate_natural_long(df: pd.DataFrame) -> pd.DataFrame:
    qc = df.copy()
    qc["dna_valid"] = qc["dna_sequence"].str.fullmatch("[ACGT]{8}")
    qc["rc_source_matches"] = (
        qc["reverse_complement_from_source"].isna()
        | qc["reverse_complement_from_source"].eq(qc["dna_sequence"].map(reverse_complement))
    )
    return qc


def consensus_by_experiment(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby(
            [
                "natural_protein_id",
                "protein_name",
                "construct_label",
                "canonical_rc",
                "publication",
                "database",
                "assay_type",
                "experimental_score_type",
            ],
            as_index=False,
        )
        .agg(
            dna_sequence=("dna_sequence", "first"),
            dna_length=("dna_length", "first"),
            experimental_score=("experimental_score", "mean"),
            median_intensity=("median_intensity", "mean"),
            z_score=("z_score", "mean"),
            n_rows_collapsed=("dna_sequence", "size"),
            experiment_ids=("experiment_id", lambda values: ";".join(sorted(set(values)))),
            replicate_ids=("replicate_id", lambda values: ";".join(sorted(set(values)))),
            source_files=("source_file", lambda values: ";".join(sorted(set(values)))),
        )
    )
    grouped["experimental_rank"] = grouped.groupby("natural_protein_id")["experimental_score"].rank(
        method="average", ascending=False
    )
    grouped["experimental_percentile"] = grouped.groupby("natural_protein_id")["experimental_score"].rank(
        pct=True, ascending=True
    )
    grouped["within_protein_normalized_rank"] = 1.0 - (
        (grouped["experimental_rank"] - 1.0)
        / grouped.groupby("natural_protein_id")["experimental_rank"].transform("max").clip(lower=1.0)
    )
    return grouped
