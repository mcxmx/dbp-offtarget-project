from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import pandas as pd
import requests

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.natural_pbm import (
    UNIPROBE_DOWNLOADS_URL,
    download_file,
    is_candidate_uniprobe_8mer_entry,
    parse_uniprobe_downloads_page,
    sha256_file,
    zip_entries,
)
from src.utils import ensure_dir, project_root


ROOT = project_root()
RAW = ensure_dir(ROOT / "data" / "raw" / "v0_4_1")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4_1")
DOCS = ensure_dir(ROOT / "docs" / "v0_4_1")
LOGS = ensure_dir(ROOT / "logs" / "v0_4_1")
RETRIEVAL_DATE = date.today().isoformat()

MIN_PROFILES = 60
MAX_ZIPS = 20
MAX_TOTAL_BYTES = 350 * 1024 * 1024
MAX_SINGLE_ZIP_BYTES = 90 * 1024 * 1024


def head_size(url: str) -> tuple[int | None, str]:
    try:
        response = requests.head(url, allow_redirects=True, timeout=30)
        response.raise_for_status()
        size = response.headers.get("Content-Length")
        return int(size) if size is not None else None, "head_ok"
    except Exception as exc:  # noqa: BLE001 - provenance audit should keep failures.
        return None, f"head_failed:{type(exc).__name__}:{exc}"


def download_uniprobe_sources() -> pd.DataFrame:
    page_path = RAW / "uniprobe_downloads_page.html"
    if not page_path.exists():
        download_file(UNIPROBE_DOWNLOADS_URL, page_path)
    html = page_path.read_text(encoding="utf-8", errors="replace")
    pubs = parse_uniprobe_downloads_page(html)
    manifest_rows = []
    total_downloaded = 0
    total_profiles = 0
    downloaded_zips = 0
    for pub in pubs:
        filename = f"{pub.dataset_code}_contig8mers.zip"
        local_path = RAW / filename
        size, head_status = head_size(pub.contig8mer_url)
        status = "not_selected"
        n_entries = 0
        sha = ""
        if local_path.exists() and local_path.stat().st_size > 0:
            try:
                status = "exists"
                sha = sha256_file(local_path)
                entries = zip_entries(local_path)
                n_entries = len([entry for entry in entries if is_candidate_uniprobe_8mer_entry(entry.filename)])
                total_profiles += n_entries
                total_downloaded += local_path.stat().st_size
            except Exception as exc:  # noqa: BLE001
                status = f"existing_file_invalid:{type(exc).__name__}:{exc}"
        elif (
            total_profiles < MIN_PROFILES
            and downloaded_zips < MAX_ZIPS
            and (size is None or size <= MAX_SINGLE_ZIP_BYTES)
            and (size is None or total_downloaded + size <= MAX_TOTAL_BYTES)
        ):
            try:
                status = download_file(pub.contig8mer_url, local_path)
                sha = sha256_file(local_path)
                entries = zip_entries(local_path)
                n_entries = len([entry for entry in entries if is_candidate_uniprobe_8mer_entry(entry.filename)])
                total_profiles += n_entries
                total_downloaded += local_path.stat().st_size
                downloaded_zips += 1
                status = "downloaded" if status == "downloaded" else "exists"
            except Exception as exc:  # noqa: BLE001
                status = f"download_failed:{type(exc).__name__}:{exc}"
        elif size is not None and size > MAX_SINGLE_ZIP_BYTES:
            status = "skipped_large_single_zip"
        elif total_downloaded + (size or 0) > MAX_TOTAL_BYTES:
            status = "skipped_total_size_cap"
        elif total_profiles >= MIN_PROFILES:
            status = "skipped_after_min_profiles_reached"
        elif downloaded_zips >= MAX_ZIPS:
            status = "skipped_after_max_zip_count"

        manifest_rows.append(
            {
                "source_database": "UniPROBE",
                "publication": pub.publication,
                "dataset_code": pub.dataset_code,
                "filename": filename,
                "url": pub.contig8mer_url,
                "local_relative_path": str(local_path.relative_to(ROOT)).replace("\\", "/") if local_path.exists() else "",
                "sha256": sha,
                "size": int(local_path.stat().st_size) if local_path.exists() else size,
                "n_8mer_profile_files": n_entries,
                "download_status": status,
                "head_status": head_status,
                "retrieval_date": RETRIEVAL_DATE,
            }
        )
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(METADATA / "natural_pbm_files.csv", index=False)
    return manifest


def write_docs(manifest: pd.DataFrame) -> None:
    downloaded = manifest[manifest["download_status"].isin(["downloaded", "exists"])]
    total_profiles = int(downloaded["n_8mer_profile_files"].sum())
    selected_codes = ", ".join(downloaded["dataset_code"].tolist())
    source_audit = f"""# v0.4.1 Natural PBM Source Audit

Audit date: {RETRIEVAL_DATE}

## Candidate Sources

| database | assay_type | number_of_proteins | protein_sequence_availability | kmer_length | score_type | replicate_availability | raw_probe_availability | species | family_diversity | downloadability | license | recommended_use |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| UniPROBE | universal PBM | publication-specific; v0.4.1 downloaded {total_profiles} profile files | factor names in profile paths; construct sequences are not generally embedded in contig8mer files | 8 | contiguous 8-mer E-score plus median/z-score where available | replicate IDs inferable for some entries | deBruijn probe files available but license-protected | mixed by publication | broad across UniPROBE publications | direct publication zip URLs work | academic research use license for files containing PBM probe sequences; E-score profiles are publicly downloadable from UniPROBE | PRIMARY SOURCE for v0.4.1 natural PBM pilot/training benchmark |
| CIS-BP 3.10 | PBM/PWM compiled specificity catalog | large database | TF information and protein metadata available in database downloads | mostly 8-mer E-scores/PWMs | E-score, Z-score, intensity, motif/PWM | varies | available in bulk archives | broad | broad | full E-score archive is ~1.19 GB; species POST archive was tested but large download was interrupted in this environment | CIS-BP database terms/citation required | SECONDARY SOURCE for future sequence metadata and larger natural PBM benchmark |
| DREAM/PBM challenge-style datasets | PBM benchmark/challenge | challenge-specific | varies | probe or k-mer level depending release | challenge-specific | varies | varies | mainly natural TFs | moderate | not selected in this pass | dataset-specific | Future validation source after harmonization audit |

## v0.4.1 Selection

PRIMARY SOURCE: UniPROBE publication-level contiguous 8-mer E-score files.

Selected UniPROBE dataset codes: {selected_codes}

Rationale: the files are directly downloadable, have explicit contiguous 8-mer E-score tables, and each profile contains reverse-complement-paired 8-mer rows suitable for per-protein ranking. This round does not claim construct-level protein sequence recovery from UniPROBE contig8mer files.

## Important Limitation

UniPROBE contig8mer profile files do not reliably provide experimental construct sequences. v0.4.1 therefore treats protein sequence recovery as a separate provenance task. Reference full-length sequences, when later added from UniProt, must be flagged as `sequence_match_to_assay=false` unless construct sequence evidence is available.
"""
    (DOCS / "NATURAL_PBM_SOURCE_AUDIT.md").write_text(source_audit, encoding="utf-8")

    harmonization = """# v0.4.1 K-mer Length Harmonization

## Natural PBM

The v0.4.1 natural PBM source uses UniPROBE contiguous 8-mer E-scores. Each source row contains an 8-mer and its reverse complement plus one E-score, so the natural sequence unit is an 8-mer reverse-complement equivalence class.

## Designed uPBM

The v0.3/v0.3.1 designed DBP benchmark from GSE237017 uses processed 7-mer PBM E-scores. Its independent DNA units are 7-mer reverse-complement classes.

## Score Compatibility

Both sources are PBM-derived enrichment/specificity scores, but their exact processing pipelines are not guaranteed identical. Scores are used for per-protein ranking, not cross-protein absolute affinity.

## Modeling Implication

Future protein-conditioned baselines must support variable-length DNA inputs or use a predeclared length harmonization strategy. v0.4.1 does not crop natural 8-mers into 7-mers and does not pad designed 7-mers into 8-mers.

## Confounding Risk

Natural held-out versus designed external performance can be confounded by k-mer length, array design, score processing, protein construct differences, and assay protocol. A drop on designed DBPs cannot by itself be interpreted as pure biological OOD.
"""
    (DOCS / "KMER_LENGTH_HARMONIZATION.md").write_text(harmonization, encoding="utf-8")

    rederive = """# Natural 7-mer Rederivation Plan

v0.4.1 does not rederive natural 7-mer scores from UniPROBE probe-level data.

Scientific reason: direct 8-mer to 7-mer truncation would merge multiple 8-mer contexts and change the PBM score definition. A valid 7-mer rederivation would need normalized probe-level data, probe sequence design, background model choices, and a documented PBM enrichment pipeline.

Future feasible plan:

1. Select one UniPROBE publication with downloadable deBruijn probe sequences and normalized probe intensity files.
2. Recompute 7-mer probe occurrence features with reverse-complement equivalence.
3. Define enrichment using a fixed robust statistic before looking at designed benchmark performance.
4. Compare rederived 8-mer scores to official UniPROBE 8-mer scores as a sanity check.
5. Only then export rederived 7-mer profiles.
"""
    (DOCS / "NATURAL_7MER_REDERIVATION_PLAN.md").write_text(rederive, encoding="utf-8")


def main() -> None:
    manifest = download_uniprobe_sources()
    write_docs(manifest)
    downloaded = manifest[manifest["download_status"].isin(["downloaded", "exists"])]
    print(downloaded[["dataset_code", "n_8mer_profile_files", "size", "download_status"]].to_string(index=False))
    print(f"Downloaded/available profile files: {int(downloaded['n_8mer_profile_files'].sum())}")


if __name__ == "__main__":
    main()
