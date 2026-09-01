from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
SOURCE_WORKBOOK = ROOT / "data" / "raw" / "v0_3_1" / "41594_2025_1669_MOESM20_ESM.xls"
METADATA_DIR = ROOT / "metadata" / "v0_3_1"
TABLES_DIR = ensure_dir(ROOT / "results" / "v0_3_1" / "tables")
DOCS_DIR = ensure_dir(ROOT / "docs" / "v0_3_1")
TOLERANCE_PERCENTILE_POINTS = 2.0


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(col) for col in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in frame.columns) + " |")
    return "\n".join(lines)


PAPER_REPORTED_PERCENTILES = {
    "DBP1": 33.19,
    "DBP3": 46.58,
    "DBP5": 86.54,
    "DBP6": 99.54,
    "DBP9": 99.89,
    "DBP35": 81.88,
    "DBP48": 97.59,
}


def protein_from_sheet(sheet_name: str) -> str:
    if "DB1" in sheet_name and "DBP1" not in sheet_name:
        return "DBP1"
    match = re.search(r"DBP(\d+)", sheet_name)
    if not match:
        raise ValueError(f"Could not parse protein ID from sheet name: {sheet_name}")
    return f"DBP{int(match.group(1))}"


def replicate_from_sheet(sheet_name: str, columns: list[str]) -> str:
    lower = sheet_name.lower()
    if "_r2" in lower or any("e-score 2" in str(col).lower() for col in columns):
        return "2"
    return "1"


def parse_source_data() -> dict[str, dict[str, pd.DataFrame]]:
    workbook = pd.ExcelFile(SOURCE_WORKBOOK)
    templates: dict[str, pd.DataFrame] = {}
    parsed: dict[str, dict[str, pd.DataFrame]] = {}
    for sheet_name in workbook.sheet_names:
        frame = pd.read_excel(SOURCE_WORKBOOK, sheet_name=sheet_name)
        columns = [str(col) for col in frame.columns]
        protein_id = protein_from_sheet(sheet_name)
        replicate = replicate_from_sheet(sheet_name, columns)
        if len(columns) >= 3:
            score_col = next(col for col in frame.columns if "E-score" in str(col))
            sub = pd.DataFrame(
                {
                    "row_index": range(len(frame)),
                    "seven_mer_a": frame.iloc[:, 0].astype(str).str.upper(),
                    "seven_mer_b": frame.iloc[:, 1].astype(str).str.upper(),
                    "e_score": pd.to_numeric(frame[score_col], errors="raise"),
                }
            )
            templates[protein_id] = sub[["row_index", "seven_mer_a", "seven_mer_b"]].copy()
        else:
            if protein_id not in templates:
                raise ValueError(f"No 7-mer template sheet found before replicate-only sheet {sheet_name}")
            score_col = frame.columns[0]
            sub = templates[protein_id].copy()
            sub["e_score"] = pd.to_numeric(frame[score_col], errors="raise")
        parsed.setdefault(protein_id, {})[replicate] = sub
    return parsed


def motif_percentile_for_replicate(frame: pd.DataFrame, motif: str) -> tuple[float, int]:
    percentile = frame["e_score"].rank(method="average", pct=True) * 100.0
    matches = frame["seven_mer_a"].str.contains(motif, regex=False) | frame["seven_mer_b"].str.contains(motif, regex=False)
    if not matches.any():
        return float("nan"), 0
    return float(percentile[matches].mean()), int(matches.sum())


def reproduce() -> pd.DataFrame:
    target_definitions = pd.read_csv(METADATA_DIR / "designed_dbp_target_definitions.csv")
    source_data = parse_source_data()
    rows = []
    for _, target_row in target_definitions.sort_values("protein_id").iterrows():
        protein_id = target_row["protein_id"]
        motif = target_row["designed_binding_site_motif"]
        replicate_values = []
        replicate_counts = []
        for replicate, frame in sorted(source_data[protein_id].items()):
            value, count = motif_percentile_for_replicate(frame, motif)
            replicate_values.append(value)
            replicate_counts.append(count)
        reproduced = float(pd.Series(replicate_values).mean())
        reported = PAPER_REPORTED_PERCENTILES[protein_id]
        diff = abs(reproduced - reported)
        rows.append(
            {
                "protein_id": protein_id,
                "paper_reported_percentile": reported,
                "our_reproduced_percentile": reproduced,
                "absolute_difference": diff,
                "motif_sequence": motif,
                "motif_length": int(target_row["designed_binding_site_motif_length"]),
                "n_matching_7mers": ";".join(str(x) for x in replicate_counts),
                "n_replicates_used": len(replicate_values),
                "replicate_percentiles": ";".join(f"{x:.6f}" for x in replicate_values),
                "matching_rule": "8192 RC-class source-data rows where either 7-mer column contains motif_sequence",
                "ranking_rule": "rank E-score within each replicate source-data sheet; higher E-score gets higher percentile; average replicate motif percentiles",
                "qc_tolerance_percentile_points": TOLERANCE_PERCENTILE_POINTS,
                "reproduction_status": "PASS" if diff <= TOLERANCE_PERCENTILE_POINTS else "FAIL",
                "notes": target_row["notes"],
            }
        )
    return pd.DataFrame(rows)


def write_notes(reproduction: pd.DataFrame) -> None:
    status = "PASS" if (reproduction["reproduction_status"] == "PASS").all() else "PARTIAL PASS"
    table = reproduction[
        [
            "protein_id",
            "paper_reported_percentile",
            "our_reproduced_percentile",
            "absolute_difference",
            "motif_sequence",
            "motif_length",
            "n_matching_7mers",
            "reproduction_status",
        ]
    ].copy()
    for col in ["paper_reported_percentile", "our_reproduced_percentile", "absolute_difference"]:
        table[col] = table[col].map(lambda value: f"{value:.3f}")
    text = f"""# Paper uPBM Reproduction Notes

Audit date: 2026-09-01

Primary source: Nature article "Computational design of sequence-specific DNA-binding proteins", DOI 10.1038/s41594-025-01669-4.

Additional sources:

- GEO GSE237017 sample metadata and processed uPBM files.
- Source Data Extended Data Fig. 8 (`41594_2025_1669_MOESM20_ESM.xls`).
- Source Data Fig. 4 (`41594_2025_1669_MOESM12_ESM.xlsx`).
- Supplementary Tables 1-3 (`41594_2025_1669_MOESM3_ESM.xlsx`).

## Paper Definition Captured

The paper reports that uPBM E-scores were used to evaluate whether 7-mers containing the designed binding-site motif were enriched among high-scoring sequences. Extended Data Fig. 8 reports motif percentile values for DBP6, DBP9, DBP48, DBP5, DBP35, DBP1, and DBP3.

The Methods/GEO processing notes describe processed PBM E-scores computed with Seed-and-wobble from Alexa 488 signal after position adjustment. These are processed experimental uPBM specificity/enrichment scores, not Kd, free energy, or binding probability.

## v0.3.1 Reproduction Rule

The Extended Data Fig. 8 source-data workbook stores 8192 rows. Each row has two reverse-complement 7-mer columns and one E-score per replicate. v0.3.1 therefore treats each source-data row as one reverse-complement equivalence class.

For each DBP:

1. Rank the 8192 source-data rows within each replicate by E-score.
2. Convert rank to percentile, with higher E-score giving higher percentile.
3. Select rows where either 7-mer column contains the designed motif.
4. Compute the mean percentile of selected rows per replicate.
5. Average replicate-level motif percentiles when replicates are available.

Tolerance: <= 2 percentile points absolute difference from the paper-reported value.

## DBP48 Target Clarification

Supplementary Table 1 records DBP48 as originally designed against target `I_b`, whose top strand is `CGCCCAAAGCCGCG`. The Fig. 4 caption states that DBP48 was analyzed with sequence C because of improved binding signal and nearly identical modeled binding sites. v0.3.1 therefore records:

- Original design target: target I.
- Experimental assay/PBM evaluation reference: target C.
- PBM motif used for percentile reproduction: `CTGACG`.

This separates original design target, assay target, and PBM motif reference.

## Reproduction Result

Overall status: {status}

{markdown_table(table)}

## Direct Paper-Derived Versus Implementation Choices

Direct from paper/source data:

- GSE237017 is the uPBM accession for the designed DBPs.
- Extended Data Fig. 8 reports the DBP-specific motif percentiles.
- The source-data workbook provides replicate E-scores for 8192 7-mer reverse-complement rows.
- Supplementary Table 1/2 provides design target IDs and exact dsDNA target sequences.
- Fig. 4 describes DBP48 analysis with sequence C.

Implementation choices:

- Use the source-data row as the RC-equivalence class.
- Match motif if either 7-mer column in that row contains the motif.
- Average replicate-level motif percentiles for DBPs with two replicates.
- Use the motif mapping recorded in `metadata/v0_3_1/designed_dbp_target_definitions.csv`.
"""
    (DOCS_DIR / "PAPER_UPBM_REPRODUCTION_NOTES.md").write_text(text, encoding="utf-8")


def main() -> None:
    reproduction = reproduce()
    reproduction.to_csv(TABLES_DIR / "paper_percentile_reproduction.csv", index=False)
    summary = pd.DataFrame(
        [
            {
                "overall_status": "PASS" if (reproduction["reproduction_status"] == "PASS").all() else "PARTIAL PASS",
                "n_proteins": int(len(reproduction)),
                "n_pass": int((reproduction["reproduction_status"] == "PASS").sum()),
                "max_absolute_difference": float(reproduction["absolute_difference"].max()),
                "tolerance_percentile_points": TOLERANCE_PERCENTILE_POINTS,
            }
        ]
    )
    summary.to_csv(TABLES_DIR / "paper_percentile_reproduction_summary.csv", index=False)
    write_notes(reproduction)
    print(reproduction.to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
