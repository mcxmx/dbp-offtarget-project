from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from file_manifest import record_file


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/v1_3_external"
PROCESSED = ROOT / "data/processed"
RESULTS = ROOT / "results/v1_3"
REPORTS = ROOT / "reports/v1.3"
FIG3 = RAW / "41594_2025_1669_MOESM16_ESM.xls"
FIG9 = RAW / "41594_2025_1669_MOESM21_ESM.xls"
LEGACY_FIG3 = ROOT / "data/raw/v1_1_competition/source_data_extended_data_fig3.xls"
DOI = "10.1038/s41594-025-01669-4"
FIG3_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM16_ESM.xls"
FIG9_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM21_ESM.xls"
PROTEINS = ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP023", "DBP035", "DBP048", "DBP056", "DBP062"]
BASES = "ACGT"


def _protein(sheet: str) -> str:
    match = re.search(r"DBP(\d+)$", sheet)
    if not match:
        raise ValueError(f"No DBP identifier in sheet name: {sheet}")
    return f"DBP{int(match.group(1)):03d}"


def _parse_landscape(path: Path, sheet: str) -> pd.DataFrame:
    table = pd.read_excel(path, sheet_name=sheet)
    expected = ["position", "original_base", "new_base", "sample", "Median PE/FITC (Normalized)"]
    if list(table.columns) != expected:
        raise ValueError(f"Unexpected columns in {path.name}/{sheet}: {list(table.columns)}")
    table = table[table["position"].astype(str) != "WT"].copy()
    table["position"] = table["position"].astype(int)
    table["original_base"] = table["original_base"].astype(str).str.upper()
    table["new_base"] = table["new_base"].astype(str).str.upper()
    table["published_mean"] = table["Median PE/FITC (Normalized)"].astype(float)
    return table


def _safe_corr(x: pd.Series, y: pd.Series, kind: str) -> float:
    a, b = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return float((spearmanr if kind == "spearman" else pearsonr)(a, b).statistic)


def _position_identity(table: pd.DataFrame, value: str) -> tuple[pd.DataFrame, np.ndarray]:
    pos = table.groupby("position", as_index=False)[value].apply(lambda x: float(np.mean(np.abs(x))))
    residual = table[value] - table.groupby("position")[value].transform("mean")
    return pos, residual.to_numpy(float)


def _source_inventory() -> pd.DataFrame:
    rows = []
    for source_figure, path in [("Extended Data Fig. 3", FIG3), ("Extended Data Fig. 9", FIG9)]:
        book = pd.ExcelFile(path)
        for sheet in book.sheet_names:
            table = pd.read_excel(path, sheet_name=sheet)
            rows.append({
                "source_figure": source_figure,
                "source_file": path.relative_to(ROOT).as_posix(),
                "sheet": sheet,
                "n_rows": len(table),
                "n_columns": len(table.columns),
                "columns": "|".join(map(str, table.columns)),
                "replicate_columns_found": int(any("rep" in str(c).lower() for c in table.columns)),
            })
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "public_source_inventory.tsv", sep="\t", index=False)
    return out


def _recover_fig3() -> tuple[pd.DataFrame, pd.DataFrame]:
    official_rows, comparison_rows = [], []
    old_book = pd.ExcelFile(LEGACY_FIG3)
    old_by_protein = {_protein(sheet): sheet for sheet in old_book.sheet_names}
    for sheet in pd.ExcelFile(FIG3).sheet_names:
        protein = _protein(sheet)
        official = _parse_landscape(FIG3, sheet)
        legacy = _parse_landscape(LEGACY_FIG3, old_by_protein[protein])
        key = ["position", "original_base", "new_base"]
        joined = official[key + ["published_mean"]].merge(
            legacy[key + ["published_mean"]], on=key, how="outer", suffixes=("_official", "_legacy"), validate="one_to_one"
        )
        joined.insert(0, "protein_id", protein)
        joined["absolute_reconstruction_error"] = (
            joined["published_mean_official"] - joined["published_mean_legacy"]
        ).abs()
        comparison_rows.append(joined)
        for _, row in official.iterrows():
            official_rows.append({
                "protein_id": protein,
                "position": int(row.position),
                "wt_base": row.original_base,
                "mutant_base": row.new_base,
                "replicate_id": "published_mean_only",
                "raw_measurement": float(row.published_mean),
                "normalized_effect": -float(row.published_mean),
                "assay_condition": "1 uM biotinylated target; 8 uM competitor; without avidity",
                "source_figure": "Extended Data Fig. 3",
                "source_file": FIG3.relative_to(ROOT).as_posix(),
                "replicate_level_available": False,
            })
    official = pd.DataFrame(official_rows).sort_values(["protein_id", "position", "mutant_base"])
    comparison = pd.concat(comparison_rows, ignore_index=True).sort_values(["protein_id", "position", "new_base"])
    official.to_parquet(PROCESSED / "v1_3_competition_mutations_official_mean.parquet", index=False)
    comparison.to_csv(RESULTS / "competition_mean_reconstruction.tsv", sep="\t", index=False)

    schema = {
        "protein_id": "string", "position": "int64", "wt_base": "string", "mutant_base": "string",
        "replicate_id": "string", "raw_measurement": "float64", "normalized_effect": "float64",
        "assay_condition": "string", "source_figure": "string", "source_file": "string",
    }
    empty = pd.DataFrame({name: pd.Series(dtype=dtype) for name, dtype in schema.items()})
    empty.to_parquet(PROCESSED / "v1_3_competition_mutations_replicates.parquet", index=False)
    return official, comparison


def _recover_dbp35opt() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sheet = "Extended_Data_Figure_9_E_DBP035"
    raw = _parse_landscape(FIG9, sheet)
    target = "GCAGATCTGCACATC"
    observed = "".join(raw.drop_duplicates("position").sort_values("position")["original_base"])
    if target[: len(observed)] != observed:
        raise ValueError(f"DBP35opt source sequence {observed} is inconsistent with target {target}")
    opt = pd.DataFrame({
        "protein_id": "DBP035opt",
        "assay": "competition",
        "target_sequence": target,
        "position": raw.position.astype(int),
        "wt_base": raw.original_base,
        "mutant_base": raw.new_base,
        "replicate": "published_summary",
        "raw_measurement": raw.published_mean.astype(float),
        "normalized_effect": -raw.published_mean.astype(float),
        "assay_condition": "20 nM biotinylated target; 160 nM competitor; without avidity",
        "source_figure": "Extended Data Fig. 9e",
        "source_file": FIG9.relative_to(ROOT).as_posix(),
        "protein_mutations": "K18V;R33N;P42Q",
        "replicate_level_available": False,
    })
    opt.to_parquet(PROCESSED / "v1_3_dbp35opt_competition_mutations.parquet", index=False)

    base = pd.read_parquet(PROCESSED / "v1_3_competition_mutations.parquet")
    base = base[base.protein_id == "DBP035"][["position", "wt_base", "mutant_base", "raw_measurement", "normalized_effect"]]
    joined = base.merge(
        opt[["position", "wt_base", "mutant_base", "raw_measurement", "normalized_effect"]],
        on=["position", "wt_base", "mutant_base"], validate="one_to_one", suffixes=("_dbp35", "_dbp35opt")
    )
    joined["experimental_effect_shift"] = joined.normalized_effect_dbp35opt - joined.normalized_effect_dbp35
    joined["dbp35_identity_residual"] = joined.normalized_effect_dbp35 - joined.groupby("position").normalized_effect_dbp35.transform("mean")
    joined["dbp35opt_identity_residual"] = joined.normalized_effect_dbp35opt - joined.groupby("position").normalized_effect_dbp35opt.transform("mean")
    joined["identity_residual_shift"] = joined.dbp35opt_identity_residual - joined.dbp35_identity_residual
    joined.to_csv(RESULTS / "dbp35opt_mutation_shift.tsv", sep="\t", index=False)

    position = joined.groupby("position", as_index=False).agg(
        dbp35_sensitivity=("normalized_effect_dbp35", lambda x: float(np.mean(np.abs(x)))),
        dbp35opt_sensitivity=("normalized_effect_dbp35opt", lambda x: float(np.mean(np.abs(x)))),
        mean_effect_shift=("experimental_effect_shift", "mean"),
        mean_abs_identity_shift=("identity_residual_shift", lambda x: float(np.mean(np.abs(x)))),
    )
    position["sensitivity_change"] = position.dbp35opt_sensitivity - position.dbp35_sensitivity
    position["continuous_change_direction"] = np.where(position.sensitivity_change > 0, "strengthened", np.where(position.sensitivity_change < 0, "weakened", "unchanged"))
    position.to_csv(RESULTS / "dbp35opt_position_shift.tsv", sep="\t", index=False)

    summary = pd.DataFrame([{
        "n_matched_mutations": len(joined),
        "n_positions": joined.position.nunique(),
        "landscape_spearman": _safe_corr(joined.normalized_effect_dbp35, joined.normalized_effect_dbp35opt, "spearman"),
        "landscape_pearson": _safe_corr(joined.normalized_effect_dbp35, joined.normalized_effect_dbp35opt, "pearson"),
        "position_sensitivity_spearman": _safe_corr(position.dbp35_sensitivity, position.dbp35opt_sensitivity, "spearman"),
        "identity_residual_spearman": _safe_corr(joined.dbp35_identity_residual, joined.dbp35opt_identity_residual, "spearman"),
        "mean_abs_experimental_shift": float(np.mean(np.abs(joined.experimental_effect_shift))),
        "mean_abs_identity_residual_shift": float(np.mean(np.abs(joined.identity_residual_shift))),
        "assay_context_status": "CONFOUNDED_DIFFERENT_ABSOLUTE_CONCENTRATIONS",
        "predicted_shift_status": "NOT_EVALUABLE_NO_DBP35OPT_MODEL_PREDICTION",
    }])
    summary.to_csv(RESULTS / "dbp35opt_shift_summary.tsv", sep="\t", index=False)
    return opt, joined, position


def _write_replicate_outputs(comparison: pd.DataFrame) -> None:
    availability = pd.DataFrame({
        "protein": PROTEINS,
        "source_figure": "Extended Data Fig. 3",
        "published_mean_available": True,
        "replicate_values_available": False,
        "status": "NOT_EVALUABLE_SOURCE_DATA_MEAN_ONLY",
    })
    availability.to_csv(RESULTS / "competition_replicate_availability.tsv", sep="\t", index=False)
    ceiling = availability[["protein", "status"]].copy()
    for col in ["global_rep_rho", "global_rep_pearson", "global_rep_mae", "global_normalized_disagreement", "position_rep_rho", "identity_rep_rho", "identity_pairwise_rep_accuracy", "identity_concordance"]:
        ceiling[col] = np.nan
    ceiling.to_csv(RESULTS / "replicate_noise_ceiling.tsv", sep="\t", index=False, na_rep="NA")

    figure2 = pd.read_csv(RESULTS / "figure2_data.tsv", sep="\t")
    pbm = figure2[figure2.method == "PBM_experimental"]
    rows = []
    metrics = {
        "GLOBAL": "global_spearman",
        "POSITION": "position_spearman",
        "IDENTITY": "within_position_residual_spearman",
        "IDENTITY_PAIRWISE": "within_position_pairwise_accuracy",
    }
    for _, row in pbm.iterrows():
        for level, metric in metrics.items():
            rows.append({"protein": row.protein_id, "comparison": "PBM_vs_competition_mean", "level": level, "metric": metric, "value": row[metric], "status": "EVALUABLE_DIFFERENT_ASSAY_ESTIMANDS"})
    for protein in PROTEINS:
        for level, metric in metrics.items():
            rows.append({"protein": protein, "comparison": "competition_R1_vs_R2", "level": level, "metric": metric, "value": np.nan, "status": "NOT_EVALUABLE_PUBLIC_SOURCE_MEAN_ONLY"})
    pd.DataFrame(rows).to_csv(RESULTS / "assay_reproducibility_matrix.tsv", sep="\t", index=False, na_rep="NA")

    max_error = comparison.absolute_reconstruction_error.max()
    median_error = comparison.absolute_reconstruction_error.median()
    report = f"""# v1.3.1 Competition Replicate Analysis

## Source audit result

The official publisher workbook for Extended Data Fig. 3 contains one `Median PE/FITC (Normalized)` column per protein and no replicate identifier or replicate-valued columns. The downloadable workbook therefore does not contain the individual replicate experiments mentioned in the article caption. `data/processed/v1_3_competition_mutations_replicates.parquet` is an intentionally empty schema-only table; no replicate values were inferred from the published mean or figure pixels.

## Published-mean parity

The newly downloaded publisher workbook was joined to the repository's prior copy by protein, position, wild-type base, and mutant base. Across {len(comparison)} mutations, maximum absolute error was {max_error:.12g} and median absolute error was {median_error:.12g}. This verifies source-file parity, not replicate reconstruction.

## Noise ceiling status

Mutation-level, position-sensitivity, and within-position identity replicate agreement are `NOT_EVALUABLE_SOURCE_DATA_MEAN_ONLY`. The requested output `results/v1_3/replicate_noise_ceiling.tsv` retains all ten proteins with explicit missing metrics. Consequently, current data cannot determine whether weak identity-level model performance reflects a model-specific failure or an identity-level experimental reproducibility ceiling.

## PBM interpretation

`results/v1_3/assay_reproducibility_matrix.tsv` places PBM-versus-competition values beside explicit unavailable competition R1-versus-R2 cells at GLOBAL, POSITION, and IDENTITY levels. PBM disagreement remains a mixture of assay estimand, biological/context, and technical differences; technical variance is not separately identifiable without the missing competition replicates.
"""
    (REPORTS / "V1_3_REPLICATE_ANALYSIS.md").write_text(report, encoding="utf-8")


def _write_reports(inventory: pd.DataFrame, comparison: pd.DataFrame, joined: pd.DataFrame, position: pd.DataFrame, hashes: dict) -> None:
    provenance = f"""# v1.3.1 External Data Provenance

Retrieval date: 2026-09-17. Paper DOI: [{DOI}](https://doi.org/{DOI}). Files were downloaded from the publisher's official Springer Nature static-content host. SHA256 was calculated once after download and cached in `artifacts/cache/file_manifest.json`; subsequent checks use path, size, and mtime first.

| figure | original filename | publisher source | bytes | SHA256 | contents used | replicate values | transformation |
|---|---|---|---:|---|---|---|---|
| Extended Data Fig. 3 | `41594_2025_1669_MOESM16_ESM.xls` | [official workbook]({FIG3_URL}) | {hashes['fig3']['size']} | `{hashes['fig3']['hash']}` | 10 DBP competition landscapes | No: one published summary column only | Mutation rows parsed; `normalized_effect = -raw_measurement`; prior mean parity checked |
| Extended Data Fig. 9 | `41594_2025_1669_MOESM21_ESM.xls` | [official workbook]({FIG9_URL}) | {hashes['fig9']['size']} | `{hashes['fig9']['hash']}` | panel 9e DBP35 K18V/R33N/P42Q landscape | No replicate-valued columns in panel 9e | Mutation rows parsed; `normalized_effect = -raw_measurement`; other raw-event sheets retained untouched |

The Extended Data Fig. 3 sheet names unexpectedly use `Extended_Data_Figure_1_*`, but their ten DBP identifiers and values match the previously archived Extended Data Fig. 3 workbook. No screenshot digitization or OCR was used. The large Fig. 9 workbook is retained verbatim; only panel 9e was converted for this benchmark.
"""
    (REPORTS / "V1_3_EXTERNAL_DATA_PROVENANCE.md").write_text(provenance, encoding="utf-8")

    summary = pd.read_csv(RESULTS / "dbp35opt_shift_summary.tsv", sep="\t").iloc[0]
    strengthened = ", ".join(map(str, position.loc[position.sensitivity_change > 0, "position"].astype(int))) or "none"
    weakened = ", ".join(map(str, position.loc[position.sensitivity_change < 0, "position"].astype(int))) or "none"
    dbp_report = f"""# v1.3.1 DBP35 to DBP35opt Analysis

## DATA-AVAILABILITY AMENDMENT

The frozen v1.3 evaluation contract pre-registered the real-perturbation endpoint but marked it unavailable. Official Extended Data Fig. 9 source data now make the experimental landscape component evaluable. This is a data-availability amendment, not a metric amendment: effect direction, `mean_abs_effect` position sensitivity, and within-position centering remain unchanged.

## Recovered data

DBP35opt is DBP35 K18V/R33N/P42Q. The official panel 9e sheet contains {len(joined)} matched substitutions across {joined.position.nunique()} positions. It provides a published summary column but no replicate-valued columns. DBP35opt used 20 nM biotinylated target and 160 nM competitor, whereas the standard DBP35 competition landscape used 1 uM target and 8 uM competitor. The comparisons below are therefore continuous descriptive comparisons with an assay-context confound, not an isolated causal protein-mutation effect.

## Experimental decomposition

- Global landscape conservation: Spearman {summary.landscape_spearman:.3f}; Pearson {summary.landscape_pearson:.3f}.
- Position-sensitivity conservation: Spearman {summary.position_sensitivity_spearman:.3f}.
- Within-position identity conservation: residual Spearman {summary.identity_residual_spearman:.3f}.
- Mean absolute mutation-landscape shift: {summary.mean_abs_experimental_shift:.3f} normalized PE/FITC units.
- Mean absolute identity-residual shift: {summary.mean_abs_identity_residual_shift:.3f}.
- Positive continuous sensitivity change at positions: {strengthened}.
- Negative continuous sensitivity change at positions: {weakened}.

No threshold was chosen and no binary gained/lost-specificity claim is made. Exact mutation and position values are in `results/v1_3/dbp35opt_mutation_shift.tsv` and `results/v1_3/dbp35opt_position_shift.tsv`.

## Model endpoint

The predicted landscape-shift endpoint remains `NOT_EVALUABLE_NO_DBP35OPT_MODEL_PREDICTION`: the recovered workbook supplies experiment, not a DBP35opt structure/model prediction. It is therefore not yet possible to say whether DeepPBS detects this real protein perturbation.
"""
    (REPORTS / "V1_3_DBP35OPT_ANALYSIS.md").write_text(dbp_report, encoding="utf-8")


def main() -> None:
    for directory in [PROCESSED, RESULTS, REPORTS, ROOT / "artifacts/cache"]:
        directory.mkdir(parents=True, exist_ok=True)
    for path in [FIG3, FIG9, LEGACY_FIG3]:
        if not path.exists():
            raise FileNotFoundError(path)
    hashes = {
        "fig3": record_file(FIG3, FIG3_URL),
        "fig9": record_file(FIG9, FIG9_URL),
    }
    inventory = _source_inventory()
    _, comparison = _recover_fig3()
    _, joined, position = _recover_dbp35opt()
    _write_replicate_outputs(comparison)
    _write_reports(inventory, comparison, joined, position, hashes)
    state = {
        "inputs": {p.relative_to(ROOT).as_posix(): {"size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns} for p in [FIG3, FIG9, LEGACY_FIG3]},
        "status": "complete",
        "official_competition_rows": 414,
        "replicate_rows_recovered": 0,
        "dbp35opt_rows": len(joined),
    }
    (RESULTS / "v1_3_1_recovery_state.json").write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(state, indent=2))


if __name__ == "__main__":
    main()
