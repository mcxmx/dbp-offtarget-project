from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.utils import ensure_dir, project_root


ROOT = project_root()
PROCESSED_V03 = ROOT / "data" / "processed" / "v0_3"
PROCESSED_DIR = ensure_dir(ROOT / "data" / "processed" / "v0_3_1")
METADATA_DIR = ROOT / "metadata" / "v0_3_1"
TABLES_DIR = ensure_dir(ROOT / "results" / "v0_3_1" / "tables")
EXPECTED_ORIENTED_7MERS = 4**7
EXPECTED_RC_CLASSES = EXPECTED_ORIENTED_7MERS // 2


def build_oriented_benchmark() -> pd.DataFrame:
    benchmark = pd.read_parquet(PROCESSED_V03 / "designed_dbp_upbm_v0_3.parquet").copy()
    target_definitions = pd.read_csv(METADATA_DIR / "designed_dbp_target_definitions.csv")
    benchmark["oriented_7mer"] = benchmark["dna_7mer"]
    benchmark["reverse_complement_7mer"] = benchmark["oriented_7mer"].map(reverse_complement)
    benchmark["canonical_7mer"] = benchmark["oriented_7mer"].map(canonical_rc)
    benchmark["is_canonical_orientation"] = benchmark["oriented_7mer"] == benchmark["canonical_7mer"]
    benchmark["experimental_escore_consensus"] = benchmark["e_score_mean"]
    benchmark["experimental_score_consensus"] = benchmark["experimental_escore_consensus"]
    benchmark["experimental_score_primary"] = benchmark["experimental_escore_consensus"]
    benchmark["experimental_score_type"] = "processed PBM E-score consensus"
    benchmark["score_processing_level"] = "GEO processed uPBM 7-mer file; replicate consensus mean within protein and oriented 7-mer"
    benchmark["experimental_unit_note"] = "oriented 7-mer row; reverse-complement equivalent rows share one RC-class experimental unit"
    if "experimental_score_raw" in benchmark.columns:
        benchmark = benchmark.drop(columns=["experimental_score_raw"])
    merge_cols = [
        "protein_id",
        "original_design_target_id",
        "original_design_target",
        "experimental_assay_target_id",
        "experimental_assay_target",
        "designed_binding_site_motif",
        "designed_binding_site_motif_length",
        "pbm_evaluation_reference",
        "confidence",
    ]
    benchmark = benchmark.merge(target_definitions[merge_cols], on="protein_id", how="left")
    return benchmark.sort_values(["protein_id", "canonical_7mer", "oriented_7mer"]).reset_index(drop=True)


def build_rc_class_benchmark(oriented: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    for (protein_id, canonical), group in oriented.groupby(["protein_id", "canonical_7mer"], sort=True):
        scores = group["experimental_escore_consensus"]
        rows.append(
            {
                "protein_id": protein_id,
                "canonical_7mer": canonical,
                "oriented_7mers": ";".join(sorted(group["oriented_7mer"].unique())),
                "n_oriented_rows": int(len(group)),
                "experimental_escore_consensus": float(scores.mean()),
                "max_oriented_escore_abs_diff": float(scores.max() - scores.min()),
                "experimental_score_type": "processed PBM E-score consensus",
                "source_gse": group["source_gse"].iloc[0],
                "source_gsms": ";".join(sorted(set(";".join(group["source_gsms"]).split(";")))),
                "n_replicates": int(group["n_replicates"].max()),
                "protein_sequence": group["protein_sequence"].iloc[0],
                "original_design_target": group["original_design_target"].iloc[0],
                "experimental_assay_target": group["experimental_assay_target"].iloc[0],
                "designed_binding_site_motif": group["designed_binding_site_motif"].iloc[0],
                "pbm_evaluation_reference": group["pbm_evaluation_reference"].iloc[0],
            }
        )
    rc_class = pd.DataFrame(rows)
    coverage_rows = []
    for protein_id, group in oriented.groupby("protein_id"):
        rc_group = rc_class[rc_class["protein_id"] == protein_id]
        coverage_rows.append(
            {
                "protein_id": protein_id,
                "n_oriented_rows": int(len(group)),
                "n_unique_oriented_7mers": int(group["oriented_7mer"].nunique()),
                "n_rc_equivalence_classes": int(group["canonical_7mer"].nunique()),
                "expected_oriented_7mers": EXPECTED_ORIENTED_7MERS,
                "expected_rc_equivalence_classes": EXPECTED_RC_CLASSES,
                "n_rc_classes_with_two_oriented_rows": int((rc_group["n_oriented_rows"] == 2).sum()),
                "n_bad_rc_classes": int((rc_group["max_oriented_escore_abs_diff"] > 1e-12).sum()),
                "max_rc_class_score_abs_diff": float(rc_group["max_oriented_escore_abs_diff"].max()),
                "coverage_status": "PASS"
                if len(group) == EXPECTED_ORIENTED_7MERS
                and group["canonical_7mer"].nunique() == EXPECTED_RC_CLASSES
                and (rc_group["max_oriented_escore_abs_diff"] <= 1e-12).all()
                else "FAIL",
            }
        )
    coverage = pd.DataFrame(coverage_rows)
    unit_summary = pd.DataFrame(
        [
            {
                "n_proteins": int(oriented["protein_id"].nunique()),
                "n_oriented_rows_total": int(len(oriented)),
                "n_oriented_7mers_per_protein": EXPECTED_ORIENTED_7MERS,
                "n_rc_equivalence_classes_total": int(len(rc_class)),
                "n_rc_equivalence_classes_per_protein": EXPECTED_RC_CLASSES,
                "independence_note": "114688 oriented rows collapse to 57344 protein-RC-class experimental units; rows are not all independent sequence units.",
            }
        ]
    )
    return rc_class, coverage, unit_summary


def main() -> None:
    oriented = build_oriented_benchmark()
    rc_class, coverage, unit_summary = build_rc_class_benchmark(oriented)
    oriented.to_parquet(PROCESSED_DIR / "designed_dbp_upbm_oriented_v0_3_1.parquet", index=False)
    oriented.head(200).to_csv(PROCESSED_DIR / "designed_dbp_upbm_oriented_v0_3_1_preview.csv", index=False)
    rc_class.to_parquet(PROCESSED_DIR / "designed_dbp_upbm_rc_class_v0_3_1.parquet", index=False)
    rc_class.head(200).to_csv(PROCESSED_DIR / "designed_dbp_upbm_rc_class_v0_3_1_preview.csv", index=False)
    coverage.to_csv(TABLES_DIR / "rc_class_coverage_qc.csv", index=False)
    unit_summary.to_csv(TABLES_DIR / "benchmark_independent_units_summary.csv", index=False)
    print(unit_summary.to_string(index=False))
    print(coverage.to_string(index=False))


if __name__ == "__main__":
    main()
