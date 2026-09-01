from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.sequence_equivalence import reverse_complement
from src.utils import ensure_dir, project_root


ROOT = project_root()
METADATA_V03 = ROOT / "metadata" / "v0_3"
METADATA_DIR = ensure_dir(ROOT / "metadata" / "v0_3_1")
RAW_DIR = ROOT / "data" / "raw"
PAPER_SOURCE_DIR = RAW_DIR / "v0_3_1"
RETRIEVAL_DATE = "2026-09-01"


PAPER_SOURCE_FILES = [
    {
        "source_type": "Nature source data",
        "source_id": "Source Data Fig. 4",
        "filename": "41594_2025_1669_MOESM12_ESM.xlsx",
        "download_url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM12_ESM.xlsx",
        "notes": "Source data for Fig. 4; used to cross-check DBP6/DBP48 uPBM and orthogonality context.",
    },
    {
        "source_type": "Nature source data",
        "source_id": "Source Data Extended Data Fig. 8",
        "filename": "41594_2025_1669_MOESM20_ESM.xls",
        "download_url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM20_ESM.xls",
        "notes": "Source data for Extended Data Fig. 8; used for published uPBM motif percentile reproduction.",
    },
]


TARGET_DEFINITIONS = [
    {
        "protein_id": "DBP1",
        "original_design_target_id": "A_c",
        "original_design_target": "TAGCAGGATGTGT",
        "experimental_assay_target_id": "A_c",
        "experimental_assay_target": "TAGCAGGATGTGT",
        "designed_binding_site_motif": "GCAGG",
        "motif_source": "Nature Extended Data Fig. 8 caption plus v0.3.1 reproduction from Source Data Extended Data Fig. 8",
        "pbm_evaluation_reference": "A-target motif GCAGG; RC-class source-data rows where either 7-mer column contains the motif",
        "confidence": "medium",
        "notes": "The paper reports low uPBM specificity percentile for DBP1. Exact motif position is not machine-readable in the source-data workbook; this motif choice reproduces the reported percentile within tolerance.",
    },
    {
        "protein_id": "DBP3",
        "original_design_target_id": "A_c",
        "original_design_target": "TAGCAGGATGTGT",
        "experimental_assay_target_id": "A_c",
        "experimental_assay_target": "TAGCAGGATGTGT",
        "designed_binding_site_motif": "GCAGGA",
        "motif_source": "Nature Extended Data Fig. 8 caption plus v0.3.1 reproduction from Source Data Extended Data Fig. 8",
        "pbm_evaluation_reference": "A-target motif GCAGGA; RC-class source-data rows where either 7-mer column contains the motif",
        "confidence": "medium",
        "notes": "The source-data workbook does not explicitly encode the motif label; this target-A motif reproduces the reported percentile within the predeclared 2 percentage-point tolerance.",
    },
    {
        "protein_id": "DBP5",
        "original_design_target_id": "B_c",
        "original_design_target": "GCAGATCTGCACATC",
        "experimental_assay_target_id": "B_c",
        "experimental_assay_target": "GCAGATCTGCACATC",
        "designed_binding_site_motif": "TGCACA",
        "motif_source": "Nature Extended Data Fig. 8 caption and target B sequence in Supplementary Table 2",
        "pbm_evaluation_reference": "B-target motif TGCACA; RC-class source-data rows where either 7-mer column contains the motif",
        "confidence": "high",
        "notes": "Motif is a target-B binding-site segment used for uPBM percentile reproduction.",
    },
    {
        "protein_id": "DBP6",
        "original_design_target_id": "B_b",
        "original_design_target": "GCAGATCTGCACATC",
        "experimental_assay_target_id": "B_b",
        "experimental_assay_target": "GCAGATCTGCACATC",
        "designed_binding_site_motif": "TGCACA",
        "motif_source": "Nature text describing specificity for the TGCACA stretch plus Extended Data Fig. 8 source-data reproduction",
        "pbm_evaluation_reference": "B-target motif TGCACA; RC-class source-data rows where either 7-mer column contains the motif",
        "confidence": "high",
        "notes": "Supplementary Table 3 records a longer DBP6 binding-site segment CTGCACAT; v0.3.1 uses the TGCACA motif because it is the motif-level PBM reference and reproduces the reported percentile within tolerance.",
    },
    {
        "protein_id": "DBP9",
        "original_design_target_id": "B_b",
        "original_design_target": "GCAGATCTGCACATC",
        "experimental_assay_target_id": "B_b",
        "experimental_assay_target": "GCAGATCTGCACATC",
        "designed_binding_site_motif": "TGCACA",
        "motif_source": "Nature Extended Data Fig. 8 caption and target B sequence in Supplementary Table 2",
        "pbm_evaluation_reference": "B-target motif TGCACA; RC-class source-data rows where either 7-mer column contains the motif",
        "confidence": "high",
        "notes": "Motif is a target-B binding-site segment used for uPBM percentile reproduction.",
    },
    {
        "protein_id": "DBP35",
        "original_design_target_id": "B_c",
        "original_design_target": "GCAGATCTGCACATC",
        "experimental_assay_target_id": "B_c",
        "experimental_assay_target": "GCAGATCTGCACATC",
        "designed_binding_site_motif": "TGCACA",
        "motif_source": "Nature Extended Data Fig. 8 caption and target B sequence in Supplementary Table 2",
        "pbm_evaluation_reference": "B-target motif TGCACA; RC-class source-data rows where either 7-mer column contains the motif",
        "confidence": "high",
        "notes": "Motif is a target-B binding-site segment used for uPBM percentile reproduction.",
    },
    {
        "protein_id": "DBP48",
        "original_design_target_id": "I_b",
        "original_design_target": "CGCCCAAAGCCGCG",
        "experimental_assay_target_id": "C",
        "experimental_assay_target": "CGACACCTGACGCG",
        "designed_binding_site_motif": "CTGACG",
        "motif_source": "Nature Fig. 4 caption states DBP48 was analyzed with sequence C; Extended Data Fig. 8 source-data reproduction supports the C-derived motif CTGACG",
        "pbm_evaluation_reference": "C-target motif CTGACG; RC-class source-data rows where either 7-mer column contains the motif",
        "confidence": "high",
        "notes": "Original design target remains I_b from Supplementary Table 1. v0.3.1 separates it from the sequence C assay/PBM evaluation reference.",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_to_root(path_text: str) -> str:
    path = Path(path_text)
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except Exception:
        return path.as_posix()


def write_portable_manifest() -> None:
    manifest = pd.read_csv(METADATA_V03 / "gse237017_file_manifest.csv")
    manifest = manifest.copy()
    manifest["local_path"] = manifest["local_path"].map(relative_to_root)
    manifest["manifest_path_type"] = "project_relative"
    manifest.to_csv(METADATA_DIR / "gse237017_file_manifest_portable.csv", index=False)


def write_paper_source_manifest() -> None:
    rows = []
    for row in PAPER_SOURCE_FILES:
        path = PAPER_SOURCE_DIR / row["filename"]
        rows.append(
            {
                **row,
                "retrieval_date": RETRIEVAL_DATE,
                "local_path": path.relative_to(ROOT).as_posix(),
                "sha256": sha256(path) if path.exists() else "",
                "file_size_bytes": path.stat().st_size if path.exists() else pd.NA,
                "download_status": "downloaded" if path.exists() else "missing",
            }
        )
    pd.DataFrame(rows).to_csv(METADATA_DIR / "paper_source_data_manifest.csv", index=False)


def write_target_definitions() -> None:
    rows = []
    for row in TARGET_DEFINITIONS:
        original = row["original_design_target"]
        assay = row["experimental_assay_target"]
        motif = row["designed_binding_site_motif"]
        rows.append(
            {
                **row,
                "original_design_target_source": "Nature Supplementary Table 1 joined to Supplementary Table 2",
                "experimental_assay_target_source": "Nature Supplementary Table 1/2; DBP48 sequence C from Fig. 4 caption and uPBM reproduction",
                "designed_binding_site_motif_length": len(motif),
                "original_design_target_length": len(original),
                "experimental_assay_target_length": len(assay),
                "motif_reverse_complement": reverse_complement(motif),
                "source_url": "https://www.nature.com/articles/s41594-025-01669-4",
                "paper_doi": "10.1038/s41594-025-01669-4",
                "retrieval_date": RETRIEVAL_DATE,
            }
        )
    pd.DataFrame(rows).sort_values("protein_id").to_csv(METADATA_DIR / "designed_dbp_target_definitions.csv", index=False)


def main() -> None:
    write_portable_manifest()
    write_paper_source_manifest()
    write_target_definitions()
    print("wrote v0.3.1 target definitions and portable manifests")


if __name__ == "__main__":
    main()
