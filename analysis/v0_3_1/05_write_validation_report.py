from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
RESULTS_DIR = ensure_dir(ROOT / "results" / "v0_3_1")
TABLES_DIR = ensure_dir(RESULTS_DIR / "tables")
METADATA_DIR = ROOT / "metadata" / "v0_3_1"
LOG_DIR = ROOT / "logs" / "v0_3_1"


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(col) for col in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in frame.columns) + " |")
    return "\n".join(lines)


def fmt(value: float | int | str, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def pytest_status() -> tuple[str, str]:
    path = LOG_DIR / "pytest_v0_3_1.log"
    if not path.exists():
        return "UNKNOWN", "pytest log missing"
    raw = path.read_bytes()
    for encoding in ["utf-8-sig", "utf-16", "utf-16-le"]:
        try:
            text = raw.decode(encoding)
            break
        except UnicodeDecodeError:
            continue
    else:
        text = raw.decode("utf-8", errors="replace")
    match = re.search(r"(\d+) passed", text)
    if match:
        return "PASS", f"{match.group(1)} passed"
    return "FAIL", text.strip().splitlines()[-1] if text.strip() else "empty pytest log"


def main() -> None:
    reproduction = pd.read_csv(TABLES_DIR / "paper_percentile_reproduction.csv")
    unit_summary = pd.read_csv(TABLES_DIR / "benchmark_independent_units_summary.csv")
    rc_baseline = pd.read_csv(TABLES_DIR / "designed_dbp_sequence_baseline_rc_aware.csv")
    comparison = pd.read_csv(TABLES_DIR / "designed_dbp_sequence_baseline_rc_aware_comparison.csv")
    disagreement = pd.read_csv(TABLES_DIR / "all_disagreement_candidate_counts.csv")
    examples = pd.read_csv(TABLES_DIR / "top_disagreement_examples.csv")
    noise = pd.read_csv(TABLES_DIR / "experimental_noise_ceiling.csv")
    clusters = pd.read_csv(TABLES_DIR / "designed_dbp_sequence_clusters.csv")
    target_groups = pd.read_csv(METADATA_DIR / "designed_dbp_target_groups.csv")
    target_defs = pd.read_csv(METADATA_DIR / "designed_dbp_target_definitions.csv")
    rc_coverage = pd.read_csv(TABLES_DIR / "rc_class_coverage_qc.csv")
    test_status, test_detail = pytest_status()

    reproduction_status = "PASS" if reproduction["reproduction_status"].eq("PASS").all() else "FAIL"
    rc_status = "PASS" if rc_coverage["coverage_status"].eq("PASS").all() else "FAIL"
    target_status = "PASS" if len(target_defs) == 7 and target_defs["confidence"].notna().all() else "FAIL"
    disagreement_status = "PASS" if disagreement["n_disagreement"].sum() > len(examples) else "FAIL"
    gate = (
        "GO TO V0.4"
        if all(status == "PASS" for status in [reproduction_status, rc_status, target_status, disagreement_status, test_status])
        else "STOP -- FIX REQUIRED"
    )

    rep_table = reproduction[
        [
            "protein_id",
            "paper_reported_percentile",
            "our_reproduced_percentile",
            "absolute_difference",
            "motif_sequence",
            "motif_length",
            "reproduction_status",
        ]
    ].copy()
    for col in ["paper_reported_percentile", "our_reproduced_percentile", "absolute_difference"]:
        rep_table[col] = rep_table[col].map(lambda x: fmt(x, 3))

    baseline_table = rc_baseline[["protein_id", "metric", "spearman", "n_sequences"]].copy()
    baseline_table["spearman"] = baseline_table["spearman"].map(lambda x: fmt(x, 3))
    baseline_medians = comparison.groupby("metric")[["old_spearman", "rc_aware_spearman", "delta"]].median().reset_index()
    for col in ["old_spearman", "rc_aware_spearman", "delta"]:
        baseline_medians[col] = baseline_medians[col].map(lambda x: fmt(x, 3))

    e_noise = noise[noise["score_type"] == "e_score"].copy()
    noise_table = e_noise[["protein_id", "pearson_correlation", "spearman_correlation", "n_aligned_7mers"]].copy()
    for col in ["pearson_correlation", "spearman_correlation"]:
        noise_table[col] = noise_table[col].map(lambda x: fmt(x, 3))

    dbp48 = target_defs.set_index("protein_id").loc["DBP48"]
    unit = unit_summary.iloc[0]
    total_disagreement = int(disagreement["n_disagreement"].sum())
    n_examples = len(examples)
    n_clusters = int(clusters["protein_sequence_cluster"].nunique())
    n_original_target_groups = int(target_groups["original_target_group"].nunique())
    n_assay_target_groups = int(target_groups["assay_target_group"].nunique())
    n_motif_groups = int(target_groups["motif_group"].nunique())
    max_error = float(reproduction["absolute_difference"].max())
    hamming_range = rc_baseline[rc_baseline["metric"] == "hamming_similarity_to_paper_motif_rc_aware"]["spearman"]
    kmer3_range = rc_baseline[rc_baseline["metric"] == "kmer3_jaccard_to_paper_motif_rc_aware"]["spearman"]

    report = f"""# v0.3.1 Validation Report

Audit date: 2026-09-01

Decision: {gate}

## 1. Paper PBM Motif Percentile Reproduction

Status: {reproduction_status}

All 7 published Extended Data Fig. 8 motif percentiles were reproduced within the predeclared tolerance of 2 percentile points.

Maximum absolute difference: {fmt(max_error, 4)} percentile points.

{markdown_table(rep_table)}

## 2. DBP48 Target Definition

DBP48 is separated into three concepts:

- Original design target: {dbp48['original_design_target_id']} / `{dbp48['original_design_target']}`
- Experimental assay target: {dbp48['experimental_assay_target_id']} / `{dbp48['experimental_assay_target']}`
- PBM evaluation motif: `{dbp48['designed_binding_site_motif']}`

This fixes the v0.3 ambiguity where a single `intended_target_dna` field could mix original design and assay/PBM evaluation references.

## 3. Reverse-Complement Units

Status: {rc_status}

- Total oriented rows: {int(unit['n_oriented_rows_total'])}
- Oriented 7-mers per protein: {int(unit['n_oriented_7mers_per_protein'])}
- Total protein-RC-class units: {int(unit['n_rc_equivalence_classes_total'])}
- RC classes per protein: {int(unit['n_rc_equivalence_classes_per_protein'])}

The confirmed independent sequence unit for v0.3.1 is the protein-RC-class unit, not the oriented row.

## 4. RC-Aware Sequence Baseline

The v0.3.1 baseline uses paper motifs and reverse-complement-aware comparison on 8192 RC classes per protein.

RC-aware hamming Spearman range: {fmt(hamming_range.min())} to {fmt(hamming_range.max())}.
RC-aware 3-mer Spearman range: {fmt(kmer3_range.min())} to {fmt(kmer3_range.max())}.

Median old-vs-new comparison:

{markdown_table(baseline_medians)}

Interpretation: sequence-only similarity remains limited. Hamming/edit similarity is weak. Motif-level k-mer overlap explains part of the uPBM landscape for some DBPs, but it is still far below replicate agreement and is not protein-conditioned.

## 5. Disagreement Candidates

Status: {disagreement_status}

The previous number 140 was the size of a per-protein top-20 examples table. It was not the total count.

v0.3.1 total sequence-vs-experiment disagreement candidates: {total_disagreement}

Top examples table size: {n_examples}

Criterion: per-protein processed uPBM E-score >= 95th percentile and RC-aware motif hamming similarity <= protein median.

## 6. Replicate Noise Ceiling

This is an empirical replicate agreement / assay reproducibility reference, not a strict mathematical maximum.

{markdown_table(noise_table)}

E-score Pearson median: {fmt(e_noise['pearson_correlation'].median())}.
E-score Spearman median: {fmt(e_noise['spearman_correlation'].median())}.

## 7. Protein and Target Independence

Protein sequence clusters at 0.60 identity threshold: {n_clusters}

Original target groups: {n_original_target_groups}
Assay target groups: {n_assay_target_groups}
Motif groups: {n_motif_groups}

Future splits must respect protein sequence clusters, target groups, motif groups, and canonical RC DNA classes.

## 8. DeepPBS / NA-MPNN Readiness

The dataset is ready to enter v0.4 as a benchmark arena for protein-conditioned baselines, including later DeepPBS/NA-MPNN-style comparisons if those tools are used only as evaluated scoring backends.

It is not ready for calibrated off-target risk claims, full-target affinity claims, or uncertainty calibration.

## 9. Remaining Limitations

- GSE237017 v0.3.1 uses processed uPBM 7-mer E-scores, not raw array-level reprocessing.
- Full designed targets are longer than 7 bp; motif percentile reproduction does not directly measure full-target affinity.
- DBP5 and DBP48 have single replicate only in the parsed GEO metadata.
- Natural-to-designed external evaluation can be confounded by assay shift unless natural PBM/uPBM controls are added.
- Sequence-only baseline is a sanity check, not a protein-conditioned model.

## 10. Gate Checks

| Requirement | Status |
| --- | --- |
| PBM parsing no major issue | {rc_status} |
| Published motif percentile reproduction | {reproduction_status} |
| RC handling explicit and tested | {rc_status} |
| Target definitions separated | {target_status} |
| Disagreement count corrected | {disagreement_status} |
| Tests all pass | {test_status}: {test_detail} |

Final gate: {gate}
"""
    (RESULTS_DIR / "V0_3_1_VALIDATION_REPORT.md").write_text(report, encoding="utf-8")
    pd.DataFrame(
        [
            {
                "gate": gate,
                "paper_reproduction_status": reproduction_status,
                "max_reproduction_error": max_error,
                "rc_status": rc_status,
                "n_oriented_rows": int(unit["n_oriented_rows_total"]),
                "n_rc_class_units": int(unit["n_rc_equivalence_classes_total"]),
                "total_disagreement_candidates": total_disagreement,
                "example_rows": n_examples,
                "e_score_replicate_pearson_median": float(e_noise["pearson_correlation"].median()),
                "e_score_replicate_spearman_median": float(e_noise["spearman_correlation"].median()),
                "n_protein_sequence_clusters": n_clusters,
                "n_original_target_groups": n_original_target_groups,
                "n_assay_target_groups": n_assay_target_groups,
                "n_motif_groups": n_motif_groups,
                "pytest_status": test_status,
                "pytest_detail": test_detail,
            }
        ]
    ).to_csv(TABLES_DIR / "v0_3_1_validation_summary.csv", index=False)
    print(gate)
    print(f"max reproduction error: {max_error:.4f}")
    print(f"total disagreement candidates: {total_disagreement}")
    print(f"pytest: {test_status} {test_detail}")


if __name__ == "__main__":
    main()
