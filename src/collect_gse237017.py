from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import requests

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.utils import ensure_dir, load_yaml, project_root


ROOT = project_root()
CONFIG = load_yaml(ROOT / "config.yaml")
V03 = CONFIG["benchmark_v0_3"]
RAW_DIR = ensure_dir(ROOT / "data" / "raw" / "gse237017")
METADATA_DIR = ensure_dir(ROOT / "metadata" / "v0_3")
LOG_DIR = ensure_dir(ROOT / "logs" / "v0_3")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def https_url(url: str) -> str:
    return url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")


def download(url: str, local_path: Path, overwrite: bool = False) -> tuple[str, int, str]:
    local_path.parent.mkdir(parents=True, exist_ok=True)
    if local_path.exists() and not overwrite:
        return "existing", local_path.stat().st_size, sha256_file(local_path)
    response = requests.get(https_url(url), stream=True, timeout=120)
    response.raise_for_status()
    with open(local_path, "wb") as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)
    return "downloaded", local_path.stat().st_size, sha256_file(local_path)


def read_family_soft(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def soft_values(block: str, key: str) -> list[str]:
    pattern = re.compile(rf"^!{re.escape(key)}(?:_\d+)? = (.*)$", re.MULTILINE)
    return [match.group(1).strip() for match in pattern.finditer(block)]


def first_soft_value(block: str, key: str) -> str:
    values = soft_values(block, key)
    return values[0] if values else ""


def split_soft(text: str) -> tuple[dict[str, Any], list[tuple[str, str]]]:
    series_block = text.split("^SAMPLE = ", 1)[0]
    series = {
        "gse_id": first_soft_value(series_block, "Series_geo_accession"),
        "title": first_soft_value(series_block, "Series_title"),
        "summary": first_soft_value(series_block, "Series_summary"),
        "overall_design": first_soft_value(series_block, "Series_overall_design"),
        "pubmed_id": first_soft_value(series_block, "Series_pubmed_id"),
        "status": first_soft_value(series_block, "Series_status"),
        "last_update_date": first_soft_value(series_block, "Series_last_update_date"),
        "sample_ids": soft_values(series_block, "Series_sample_id"),
    }
    samples: list[tuple[str, str]] = []
    for sample_part in text.split("^SAMPLE = ")[1:]:
        gsm_id, block = sample_part.split("\n", 1)
        samples.append((gsm_id.strip(), block))
    return series, samples


def normalize_protein_id(text: str) -> str:
    match = re.search(r"\bDBP\s*(\d+)\b", text, flags=re.IGNORECASE)
    return f"DBP{match.group(1)}" if match else ""


def parse_concentration(text: str) -> str:
    match = re.search(r"(\d+(?:\.\d+)?)\s*uM", text, flags=re.IGNORECASE)
    return f"{match.group(1)}uM" if match else ""


def parse_replicate(text: str) -> str:
    match = re.search(r"replicate\s*(\d+)", text, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    match = re.search(r"-r(\d+)", text, flags=re.IGNORECASE)
    return match.group(1) if match else ""


def file_type(filename: str) -> str:
    lowered = filename.lower()
    if "_7mers_" in lowered:
        return "processed_7mer"
    if "_rawdata" in lowered:
        return "raw_spot_data"
    return "other"


def parse_samples(series: dict[str, Any], sample_blocks: list[tuple[str, str]], metadata_only: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_rows = []
    manifest_rows = []
    retrieval_date = date.today().isoformat()
    for gsm_id, block in sample_blocks:
        title = first_soft_value(block, "Sample_title")
        platform_id = first_soft_value(block, "Sample_platform_id")
        source_name = first_soft_value(block, "Sample_source_name_ch1")
        characteristics = soft_values(block, "Sample_characteristics_ch1")
        processing = soft_values(block, "Sample_data_processing")
        supplementary = soft_values(block, "Sample_supplementary_file")
        joined = " ".join([title, *characteristics, *supplementary])
        protein_id = normalize_protein_id(joined)
        protein_name = protein_id
        protein_concentration = parse_concentration(joined)
        replicate = parse_replicate(joined)
        notes = " | ".join(processing)
        sample_rows.append(
            {
                "gse_id": V03["gse_id"],
                "gsm_id": gsm_id,
                "protein_id": protein_id,
                "protein_name": protein_name,
                "protein_concentration": protein_concentration,
                "replicate": replicate,
                "sample_title": title,
                "platform_id": platform_id,
                "experiment_type": "universal protein-binding microarray (uPBM)",
                "source_name": source_name,
                "supplementary_files": ";".join(supplementary),
                "retrieval_date": retrieval_date,
                "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_id}",
                "notes": notes,
            }
        )
        for url in supplementary:
            filename = Path(url).name
            local_path = RAW_DIR / filename
            if metadata_only:
                if local_path.exists():
                    status = "existing"
                    size = local_path.stat().st_size
                    digest = sha256_file(local_path)
                else:
                    status = "metadata_only_pending_download"
                    size = pd.NA
                    digest = ""
            else:
                try:
                    status, size, digest = download(url, local_path)
                except Exception as exc:
                    status = f"error: {exc!r}"
                    size = pd.NA
                    digest = ""
            manifest_rows.append(
                {
                    "gsm_id": gsm_id,
                    "filename": filename,
                    "file_type": file_type(filename),
                    "download_url": https_url(url),
                    "local_path": str(local_path),
                    "sha256": digest,
                    "file_size_bytes": size,
                    "download_status": status,
                }
            )
    samples = pd.DataFrame(sample_rows).sort_values(["protein_id", "replicate", "gsm_id"])
    manifest = pd.DataFrame(manifest_rows).sort_values(["gsm_id", "filename"])
    return samples, manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Collect GSE237017 GEO metadata and supplementary files.")
    parser.add_argument("--metadata-only", action="store_true", help="Parse GEO metadata without downloading GSM supplementary files.")
    args = parser.parse_args(argv)

    soft_path = RAW_DIR / "GSE237017_family.soft.gz"
    filelist_path = RAW_DIR / "filelist.txt"
    soft_status, soft_size, soft_hash = download(V03["geo_family_soft_url"], soft_path)
    filelist_status, filelist_size, filelist_hash = download(V03["geo_filelist_url"], filelist_path)
    text = read_family_soft(soft_path)
    series, sample_blocks = split_soft(text)
    samples, manifest = parse_samples(series, sample_blocks, metadata_only=args.metadata_only)

    series_out = {
        **series,
        "geo_series_url": V03["geo_series_url"],
        "family_soft_url": V03["geo_family_soft_url"],
        "family_soft_local_path": str(soft_path),
        "family_soft_sha256": soft_hash,
        "family_soft_size_bytes": soft_size,
        "family_soft_download_status": soft_status,
        "filelist_url": V03["geo_filelist_url"],
        "filelist_local_path": str(filelist_path),
        "filelist_sha256": filelist_hash,
        "filelist_size_bytes": filelist_size,
        "filelist_download_status": filelist_status,
        "retrieval_date": date.today().isoformat(),
    }
    with open(METADATA_DIR / "gse237017_series_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(series_out, handle, indent=2, sort_keys=True)
    samples.to_csv(METADATA_DIR / "gse237017_samples.csv", index=False)
    manifest.to_csv(METADATA_DIR / "gse237017_file_manifest.csv", index=False)
    summary = {
        "n_samples": int(len(samples)),
        "n_proteins": int(samples["protein_id"].nunique()),
        "n_manifest_files": int(len(manifest)),
        "n_processed_7mer_files": int((manifest["file_type"] == "processed_7mer").sum()),
        "n_raw_spot_files": int((manifest["file_type"] == "raw_spot_data").sum()),
        "metadata_only": bool(args.metadata_only),
    }
    with open(METADATA_DIR / "gse237017_collection_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(samples[["gsm_id", "protein_id", "protein_concentration", "replicate", "sample_title"]].to_string(index=False))


if __name__ == "__main__":
    main()
