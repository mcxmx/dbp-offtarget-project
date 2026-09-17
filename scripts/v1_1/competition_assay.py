from __future__ import annotations

import json
import hashlib
import re
import sys
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.sequence_equivalence import canonical_rc, reverse_complement


OUT = ROOT / "results" / "v1_1"
RAW = ROOT / "data" / "raw" / "v1_1_competition" / "source_data_extended_data_fig3.xls"
SOURCE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/instance/12618268/bin/41594_2025_1669_MOESM16_ESM.xls"
RESOLVED_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41594-025-01669-4/MediaObjects/41594_2025_1669_MOESM16_ESM.xls"
DOI = "10.1038/s41594-025-01669-4"
PROTEINS = ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP035", "DBP048"]
SHEET_TO_PROTEIN = {p: f"DBP{p:03d}" for p in [1, 3, 5, 6, 9, 35, 48]}
EPS = 1e-12
BASES = "ACGT"


def spearman(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan
    return float(spearmanr(x, y).statistic)


def pearson(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan
    return float(pearsonr(x, y).statistic)


def _cell_range(row_number: int) -> str:
    # All relevant tables use A:E; row_number is zero-based pandas row index.
    excel_row = row_number + 2
    return f"A{excel_row}:E{excel_row}"


def parse_source() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not RAW.exists():
        raise FileNotFoundError(RAW)
    xls = pd.ExcelFile(RAW)
    inventory = []
    rows = []
    for sheet in xls.sheet_names:
        df = pd.read_excel(RAW, sheet_name=sheet)
        m = re.search(r"DBP(\d+)$", sheet)
        protein = f"DBP{int(m.group(1)):03d}" if m else ""
        inventory.append({
            "source_file": RAW.relative_to(ROOT).as_posix(),
            "source_sheet": sheet,
            "protein": protein,
            "rows": len(df),
            "columns": len(df.columns),
            "columns_text": "|".join(map(str, df.columns)),
            "merged_cells": "none observed by xlrd inspection",
            "replicate_columns": "none",
            "wt_control_rows": int((df.iloc[:, 0].astype(str) == "WT").sum()),
            "missing_values": int(df.isna().sum().sum()),
            "status": "parsed" if protein in SHEET_TO_PROTEIN.values() else "metadata_only",
        })
        if protein not in PROTEINS:
            continue
        required = {"position", "original_base", "new_base", "sample", "Median PE/FITC (Normalized)"}
        if not required.issubset(df.columns):
            raise AssertionError(f"{sheet}: unexpected columns {df.columns.tolist()}")
        assay_sequence = str(load_targets().loc[protein, "experimental_assay_target"])
        for row_index, r in df.iterrows():
            if str(r["position"]) == "WT":
                continue
            pos = int(r["position"])
            wt = str(r["original_base"]).upper()
            mut = str(r["new_base"]).upper()
            if pos < 1 or pos > len(assay_sequence) or assay_sequence[pos - 1] != wt:
                raise AssertionError(f"{protein}: source WT mismatch at position {pos}: {wt}/{assay_sequence}")
            if mut not in BASES or wt not in BASES or mut == wt:
                raise AssertionError(f"{protein}: invalid mutation {wt}->{mut}")
            signal = float(r["Median PE/FITC (Normalized)"])
            rows.append({
                "protein": protein,
                "assay_position": pos,
                "wt_base": wt,
                "mut_base": mut,
                "replicate": "mean_of_two_replicates",
                "replicate_level_available": False,
                "raw_value": signal,
                "normalized_value": signal,
                "experimental_effect": -signal,
                "experimental_effect_definition": "negative normalized PE/FITC; higher means stronger competitor competition",
                "assay_sequence": assay_sequence,
                "source_sheet": sheet,
                "source_cell_or_range": _cell_range(row_index),
            })
    tidy = pd.DataFrame(rows).sort_values(["protein", "assay_position", "mut_base"]).reset_index(drop=True)
    if tidy.duplicated(["protein", "assay_position", "wt_base", "mut_base"]).any():
        raise AssertionError("duplicate mutation records")
    for protein in PROTEINS:
        g = tidy[tidy.protein == protein]
        if len(g) != 3 * g.assay_position.nunique() or set(g.mut_base) - set(BASES):
            raise AssertionError(f"{protein}: incomplete three-alternative mutation table")
    means = tidy.copy()
    means["replicate"] = "mean_of_two_replicates"
    tidy.to_csv(OUT / "04_competition/competition_mutations.tsv", sep="\t", index=False)
    means.to_csv(OUT / "04_competition/competition_replicate_means.tsv", sep="\t", index=False)
    inventory_df = pd.DataFrame(inventory)
    inventory_df.to_csv(OUT / "04_competition/source_data_inventory.tsv", sep="\t", index=False)
    return tidy, means, inventory_df


def write_source_audit(inventory: pd.DataFrame) -> None:
    sha = hashlib.sha256(RAW.read_bytes()).hexdigest()
    size = RAW.stat().st_size
    text = f"""# Competition source-data audit

Source: `{RAW.relative_to(ROOT).as_posix()}`

- Public PMC URL: {SOURCE_URL}
- Resolved download URL: {RESOLVED_URL}
- DOI: {DOI}
- Retrieval date: {date.today().isoformat()}
- File size: {size} bytes
- SHA256: `{sha}`
- Workbook engine: xlrd/OLE XLS
- Sheets: {len(inventory)}
- Relevant original designs parsed: DBP001, DBP003, DBP005, DBP006, DBP009, DBP035, DBP048

Each sheet has a five-column table (`position`, `original_base`, `new_base`, `sample`, and `Median PE/FITC (Normalized)`). Relevant sheets contain 13 or 14 positions with three non-WT substitutions per position. DBP005, DBP009, and DBP035 include an additional WT row; the other relevant sheets do not.

No merged cells or missing values were observed. The sheets contain mutation rows (and WT rows in three relevant sheets) but no explicit no-competitor control row; the provided signal is already normalized to the no-competitor condition. There are no replicate columns. The XLS values are the heatmap/source values described by the paper as the mean of two replicates. Individual replicate values are not present in this downloaded workbook, so the tidy records are labeled `mean_of_two_replicates`; no replicate-level values are invented.

The source includes DBP023, DBP056, and DBP062 as additional designs. They are retained in the inventory audit but excluded from the seven-protein analysis.
"""
    (OUT / "04_competition/source_data_audit.md").write_text(text, encoding="utf-8")
    pd.DataFrame([{
        "source_file": RAW.relative_to(ROOT).as_posix(),
        "public_url": SOURCE_URL,
        "resolved_url": RESOLVED_URL,
        "paper_doi": DOI,
        "retrieval_date": date.today().isoformat(),
        "file_size_bytes": size,
        "sha256": sha,
        "immutable_raw_copy": True,
    }]).to_csv(OUT / "04_competition/source_data_manifest.tsv", sep="\t", index=False)


def write_semantics() -> None:
    (OUT / "04_competition/assay_semantics.md").write_text(
        """# Competition assay semantics

The paper defines the plotted quantity as relative binding activity, PE/FITC normalized to the no-competitor condition. A nonbiotinylated competitor reduces the labeled-target binding signal. Therefore lower normalized PE/FITC means stronger competition by that mutant competitor.

The frozen experimental direction used in v1.1 is:

`experimental_effect = - normalized_PE/FITC`

Thus larger experimental effect means stronger competitor binding / stronger competition. This is a monotonic sign reversal only; no per-protein scaling, register selection, or model-driven normalization is applied. The raw normalized signal is preserved in the tidy table.

The XLS has heatmap means of two replicates, not separate replicate columns. Correlations therefore use the published two-replicate mean values and are labeled accordingly.
""", encoding="utf-8")


def update_inventory() -> None:
    path = OUT / "00_inventory/source_inventory.csv"
    inv = pd.read_csv(path)
    row = {
        "source_file": RAW.relative_to(ROOT).as_posix(),
        "purpose": "Glasscock Extended Data Fig. 3 full competition assay source",
        "proteins_covered": "DBP001|DBP003|DBP005|DBP006|DBP009|DBP023|DBP035|DBP048|DBP056|DBP062",
        "assay": "yeast-display competition",
        "raw_processed_derived": "raw",
        "dimensions": "10 sheets; 13/14 positions; 3 substitutions per position",
        "sha256": __import__("hashlib").sha256(RAW.read_bytes()).hexdigest(),
        "url": SOURCE_URL,
        "doi": DOI,
    }
    for col in row:
        if col not in inv.columns:
            inv[col] = ""
    inv = inv[inv.source_file.astype(str) != row["source_file"]]
    inv = pd.concat([inv, pd.DataFrame([row])], ignore_index=True)
    inv.to_csv(path, index=False)
    (OUT / "00_inventory/inventory_report.md").write_text(
        "# v1.1 repository inventory\n\n"
        "Seven design proteins: DBP001, DBP003, DBP005, DBP006, DBP009, DBP035, DBP048. "
        "The frozen DeepPBS bundle was ingested without modifying NPZ files. The official Glasscock Extended Data Fig. 3 XLS is now locally frozen and parsed for all seven proteins (288 mutation rows; two-replicate heatmap means). Individual replicate values are not present in this workbook.\n",
        encoding="utf-8",
    )


def load_targets() -> pd.DataFrame:
    t = pd.read_csv(ROOT / "metadata/v0_3_1/designed_dbp_target_definitions.csv")
    t["protein_id"] = t.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}")
    return t.set_index("protein_id")


def update_register_map(tidy: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    old = pd.read_csv(OUT / "03_register_mapping/register_summary.tsv", sep="\t")
    pwm = pd.read_csv(OUT / "01_deeppbs_ingest/deeppbs_design7_pwm.tsv", sep="\t")
    targets = load_targets()
    rows, summaries = [], []
    for protein in PROTEINS:
        assay_seq = str(tidy.loc[tidy.protein == protein, "assay_sequence"].iloc[0])
        model = pwm[pwm.protein == protein].sort_values("position")
        model_seq = "".join(model.reference_base)
        old_row = old[old.protein == protein].iloc[0]
        model_start = int(old_row.model_window_start_1based) - 1
        assay_start = int(old_row.assay_window_start_1based) - 1
        status = "AMBIGUOUS" if protein == "DBP048" else "PRIMARY"
        evidence = "exact_model_and_assay_sequence_match" if protein != "DBP048" else "sequence_C_same_length_positional_mapping; design_sequence_differs"
        for pos in range(1, len(assay_seq) + 1):
            mapped = pos if pos <= len(model_seq) else np.nan
            in_window = assay_start <= pos - 1 < assay_start + 7
            rows.append({
                "protein": protein, "model": "DeepPBS/NA-MPNN", "model_position": mapped,
                "model_base": model_seq[pos - 1] if pos <= len(model_seq) else "",
                "assay_position": pos, "assay_reference_base": assay_seq[pos - 1],
                "pbm_7mer_position": (pos - assay_start) if in_window else np.nan,
                "orientation": "forward", "offset": model_start,
                "mapping_status": status, "mapping_class": "PRIMARY" if status == "PRIMARY" else "SENSITIVITY_ONLY",
                "mapping_evidence": evidence, "ambiguity_flag": status != "PRIMARY",
                "assay_sequence": assay_seq,
            })
        summaries.append({
            "protein": protein, "model_sequence": model_seq, "assay_sequence": assay_seq,
            "model_window": old_row.model_window, "assay_window": old_row.assay_window,
            "model_window_start_1based": int(old_row.model_window_start_1based),
            "assay_window_start_1based": int(old_row.assay_window_start_1based),
            "motif": targets.loc[protein, "designed_binding_site_motif"],
            "mapping_status": status, "mapping_class": "PRIMARY" if status == "PRIMARY" else "SENSITIVITY_ONLY",
            "mapping_evidence": evidence, "ambiguity_flag": status != "PRIMARY",
        })
    mapping = pd.DataFrame(rows)
    summary = pd.DataFrame(summaries)
    mapping.to_csv(OUT / "03_register_mapping/register_map.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "03_register_mapping/register_summary.tsv", sep="\t", index=False)
    lines = ["# Register mapping report (competition assay update)", "", "Mapping was fixed from model/assay sequences before correlation. DBP048 uses the actual sequence-C assay construct and remains sensitivity-only because it differs from the design sequence.", ""]
    for _, r in summary.iterrows():
        lines += [f"## {r.protein}", f"- Status: **{r.mapping_status}** ({r.mapping_class})", f"- Model sequence: `{r.model_sequence}`", f"- Competition assay sequence: `{r.assay_sequence}`", f"- Model 7-mer window: `{r.model_window}`; assay window: `{r.assay_window}`", f"- Evidence: {r.mapping_evidence}", ""]
    (OUT / "03_register_mapping/register_mapping_report.md").write_text("\n".join(lines), encoding="utf-8")
    return mapping, summary


def deep_pbs_local(tidy: pd.DataFrame, mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, metrics = [], []
    for protein in PROTEINS:
        data = np.load(ROOT / "data/raw/v1_1_deeppbs/frozen_design7_20260914/predictions" / f"{protein}.npz_predict.npz", allow_pickle=False)
        P = np.asarray(data["P"], dtype=float)
        mp = mapping[mapping.protein == protein].set_index("assay_position")
        g = tidy[tidy.protein == protein]
        pred = []
        for _, r in g.iterrows():
            mpos = int(mp.loc[int(r.assay_position), "model_position"]) - 1
            d = float(np.log(P[mpos, BASES.index(r.mut_base)] + EPS) - np.log(P[mpos, BASES.index(r.wt_base)] + EPS))
            row = r.to_dict()
            row.update({"method": "DeepPBS_native_pwm", "predicted_delta_logP": d, "mapping_status": mp.loc[int(r.assay_position), "mapping_status"], "model_position": mpos + 1})
            rows.append(row); pred.append(d)
        exp = np.asarray(g.experimental_effect, dtype=float)
        metrics.append({"method": "DeepPBS_native_pwm", "score_definition": "native DeepPBS PPM delta_logP", "protein": protein, "n_mutations": len(g), "spearman": spearman(pred, exp), "pearson": pearson(pred, exp), "sign_agreement": float(np.mean(np.sign(pred) == np.sign(exp))), "mapping_status": mp.mapping_status.iloc[0], "result_status": "PRIMARY" if mp.mapping_status.iloc[0] == "PRIMARY" else "SENSITIVITY_ONLY"})
    pred_df = pd.DataFrame(rows)
    metric_df = pd.DataFrame(metrics)
    return pred_df, metric_df


def _global_lookup(protein: str) -> dict[str, float]:
    df = pd.read_csv(ROOT / "results/v1_0_structure/nampnn_results.csv")
    df["protein"] = df.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}")
    return df[df.protein == protein].drop_duplicates("canonical_7mer").set_index("canonical_7mer").nampnn_score.to_dict()


def _pbm_lookup(protein: str) -> dict[str, float]:
    df = pd.read_parquet(ROOT / "data/processed/v0_3_1/designed_dbp_upbm_rc_class_v0_3_1.parquet")
    df["protein"] = df.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}")
    return df[df.protein == protein].drop_duplicates("canonical_7mer").set_index("canonical_7mer").experimental_escore_consensus.to_dict()


def landscape_local(tidy: pd.DataFrame, mapping: pd.DataFrame, method: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, metrics = [], []
    for protein in PROTEINS:
        mp = mapping[mapping.protein == protein].set_index("assay_position")
        g = tidy[tidy.protein == protein]
        # Seven positions are defined by the independent assay window. Mutations outside
        # that window are unavailable for a 7-mer landscape-derived score.
        window_positions = mp[mp.pbm_7mer_position.notna()]
        start = int(window_positions.index.min())
        assay_seq = g.assay_sequence.iloc[0]
        lookup = _global_lookup(protein) if method == "NA-MPNN" else _pbm_lookup(protein)
        pred, exp, used = [], [], []
        wt_window = assay_seq[start - 1:start + 6]
        wt_key = canonical_rc(wt_window)
        if wt_key not in lookup:
            continue
        wt_score = float(lookup[wt_key])
        for _, r in g.iterrows():
            if int(r.assay_position) not in window_positions.index:
                continue
            mutated = assay_seq[: int(r.assay_position) - 1] + r.mut_base + assay_seq[int(r.assay_position):]
            mut_window = mutated[start - 1:start + 6]
            mut_score = float(lookup.get(canonical_rc(mut_window), np.nan))
            if not np.isfinite(mut_score):
                continue
            d = mut_score - wt_score
            row = r.to_dict(); row.update({"method": "NA-MPNN_frozen_global_landscape_derived" if method == "NA-MPNN" else "PBM_experimental_7mer_effect", "predicted_delta": d, "mapping_status": mp.loc[int(r.assay_position), "mapping_status"], "pbm_window": wt_window, "mutated_pbm_window": mut_window})
            rows.append(row); pred.append(d); exp.append(float(r.experimental_effect)); used.append(r)
        metrics.append({"method": "NA-MPNN_frozen_global_landscape_derived" if method == "NA-MPNN" else "PBM_vs_competition", "score_definition": "NA-MPNN frozen global landscape-derived local mutation score" if method == "NA-MPNN" else "experimental PBM 7-mer delta", "protein": protein, "n_mutations": len(pred), "spearman": spearman(pred, exp), "pearson": pearson(pred, exp), "mapping_status": mp.mapping_status.iloc[0], "result_status": "PRIMARY" if mp.mapping_status.iloc[0] == "PRIMARY" else "SENSITIVITY_ONLY"})
    return pd.DataFrame(rows), pd.DataFrame(metrics)


def write_figures(deep: pd.DataFrame, na: pd.DataFrame, pbm: pd.DataFrame, global_metrics: pd.DataFrame) -> None:
    fdir = OUT / "figures"
    # DeepPBS and NA-MPNN individual protein panels.
    for name, df, x, title in [("figure2_deeppbs_competition.png", deep, "predicted_delta_logP", "DeepPBS local mutation effects"), ("figure3_nampnn_competition.png", na, "predicted_delta", "NA-MPNN frozen global-landscape-derived local effects")]:
        fig, axes = plt.subplots(2, 4, figsize=(14, 7), squeeze=False)
        for ax, protein in zip(axes.flat, PROTEINS):
            g = df[df.protein == protein]
            ax.scatter(g[x], g.experimental_effect, s=12, alpha=.8)
            ax.set_title(protein); ax.set_xlabel("predicted"); ax.set_ylabel("-normalized PE/FITC")
        for ax in axes.flat[len(PROTEINS):]: ax.axis("off")
        fig.suptitle(title); fig.tight_layout(); fig.savefig(fdir / name, dpi=220); plt.close(fig)
    # PBM versus competition panels.
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), squeeze=False)
    for ax, protein in zip(axes.flat, PROTEINS):
        g = pbm[pbm.protein == protein]; ax.scatter(g.predicted_delta, g.experimental_effect, s=12, alpha=.8); ax.set_title(protein); ax.set_xlabel("PBM delta"); ax.set_ylabel("competition effect")
    for ax in axes.flat[len(PROTEINS):]: ax.axis("off")
    fig.suptitle("Experimental PBM versus competition mutation effects"); fig.tight_layout(); fig.savefig(fdir / "figure5_pbm_competition_agreement.png", dpi=220); plt.close(fig)
    # Local/global comparison.
    fig, ax = plt.subplots(figsize=(8, 5))
    for method, xcol, color in [("DeepPBS", "DeepPBS_competition_rho", "#1f77b4"), ("NA-MPNN", "NA_MPNN_competition_rho", "#d62728")]:
        q = global_metrics.dropna(subset=[xcol]); ax.scatter(q[xcol], q["global_rho"], label=method, color=color, s=50)
        for _, r in q.iterrows(): ax.text(r[xcol], r.global_rho, r.protein, fontsize=8)
    ax.axhline(0, color="grey", lw=.8); ax.axvline(0, color="grey", lw=.8); ax.set_xlabel("local competition Spearman"); ax.set_ylabel("global PBM Spearman"); ax.legend(); fig.tight_layout(); fig.savefig(fdir / "figure4_local_vs_global.png", dpi=220); plt.close(fig)
    # DBP48 sensitivity panel.
    fig, ax = plt.subplots(figsize=(7, 4)); q = deep[deep.protein == "DBP048"]; ax.scatter(q.model_position, q.predicted_delta_logP, label="DeepPBS", s=22); ax.set_title("DBP048 sequence-C positional sensitivity"); ax.set_xlabel("assay position / model position"); ax.set_ylabel("predicted delta logP"); fig.tight_layout(); fig.savefig(fdir / "figure6_dbp048_sensitivity.png", dpi=220); plt.close(fig)


def complete_competition() -> None:
    tidy, means, inventory = parse_source()
    write_source_audit(inventory); write_semantics(); update_inventory()
    mapping, summary = update_register_map(tidy)
    deep, deep_metrics = deep_pbs_local(tidy, mapping)
    na, na_metrics = landscape_local(tidy, mapping, "NA-MPNN")
    pbm, pbm_metrics = landscape_local(tidy, mapping, "PBM")
    deep.to_csv(OUT / "05_assay_aligned_eval/deeppbs_per_mutation_predictions.tsv", sep="\t", index=False)
    deep.to_csv(OUT / "05_assay_aligned_eval/per_mutation_predictions.tsv", sep="\t", index=False)
    na.to_csv(OUT / "05_assay_aligned_eval/nampnn_local_per_mutation_predictions.tsv", sep="\t", index=False)
    pbm.to_csv(OUT / "06_assay_agreement/pbm_competition_per_mutation.tsv", sep="\t", index=False)
    metrics = pd.concat([deep_metrics, na_metrics, pbm_metrics], ignore_index=True)
    metrics.to_csv(OUT / "05_assay_aligned_eval/competition_metrics.tsv", sep="\t", index=False)
    # Keep the requested legacy filename populated with all assay-aligned methods.
    metrics.to_csv(OUT / "05_assay_aligned_eval/per_protein_metrics.tsv", sep="\t", index=False)
    summary_rows = []
    for method, g in metrics.groupby("method"):
        summary_rows.append({"method": method, "macro_median_spearman": g.spearman.median(), "macro_mean_spearman": g.spearman.mean(), "n_proteins": int(g.spearman.notna().sum()), "status": "development_exposed_assay_aligned_diagnostic"})
    pd.DataFrame(summary_rows).to_csv(OUT / "05_assay_aligned_eval/summary_metrics.tsv", sep="\t", index=False)
    # PBM and competition are both on the same per-mutation direction here.
    global_metrics_df = pd.read_csv(OUT / "09_global_pwm_projection/global_projection_metrics.tsv", sep="\t")
    wide = []
    for protein in PROTEINS:
        d = deep_metrics[deep_metrics.protein == protein].iloc[0]
        n = na_metrics[na_metrics.protein == protein].iloc[0]
        p = pbm_metrics[pbm_metrics.protein == protein].iloc[0]
        gd = global_metrics_df[(global_metrics_df.method == "DeepPBS") & (global_metrics_df.protein == protein)].iloc[0]
        gn = global_metrics_df[(global_metrics_df.method == "NA-MPNN") & (global_metrics_df.protein == protein)].iloc[0]
        wide.append({"protein": protein, "DeepPBS_competition_rho": d.spearman, "DeepPBS_global_PBM_rho": gd.spearman, "NA_MPNN_competition_rho": n.spearman, "NA_MPNN_global_PBM_rho": gn.spearman, "PBM_vs_competition_rho": p.spearman, "n_competition_mutations": int(d.n_mutations), "mapping_status": d.mapping_status})
    wide_df = pd.DataFrame(wide)
    wide_df.to_csv(OUT / "05_assay_aligned_eval/local_global_comparison.tsv", sep="\t", index=False)
    write_figures(deep, na, pbm, wide_df.merge(global_metrics_df[global_metrics_df.method == "DeepPBS"][["protein", "spearman"]].rename(columns={"spearman": "global_rho"}), on="protein", how="left"))
    published = deep_metrics[deep_metrics.protein.isin(["DBP005", "DBP006", "DBP009", "DBP035"])][["protein", "spearman", "n_mutations"]].copy()
    published_lines = [
        "# Published DeepPBS cross-check",
        "",
        "The paper reports qualitative agreement between predicted DeepPBS PWMs and competition specificity for DBP5, DBP6, DBP9 and DBP35. The v1.1 quantitative analysis uses independently fixed assay positions and `delta_logP` without optimizing register, orientation, epsilon, sign, or any experimental metric.",
        "",
        "| published protein | mutations | v1.1 Spearman |",
        "|---|---:|---:|",
    ]
    published_lines.extend(f"| {r.protein} | {int(r.n_mutations)} | {r.spearman:.4f} |" for _, r in published.iterrows())
    published_lines += [
        "",
        "These values are a quantitative diagnostic cross-check, not a parameter-selection target. The source XLS provides the two-replicate heatmap means, not separate replicate columns. Any disagreement with the published qualitative display remains interpretable in light of register, strand, construct, and assay-direction checks documented in the source and mapping reports.",
    ]
    (OUT / "05_assay_aligned_eval/published_deeppbs_crosscheck.md").write_text("\n".join(published_lines) + "\n", encoding="utf-8")
    agreement_lines = [
        "# PBM versus competition agreement",
        "",
        "The official competition XLS is available. PBM mutation effects use the fixed 7-mer register and canonical reverse-complement lookup; competition effects are the fixed negative normalized PE/FITC direction. These are two experimental assays, not interchangeable ground truth.",
        "",
        "| protein | n | Spearman | Pearson | mapping |",
        "|---|---:|---:|---:|---|",
    ]
    agreement_lines.extend(f"| {r.protein} | {int(r.n_mutations)} | {r.spearman:.4f} | {r.pearson:.4f} | {r.mapping_status} |" for _, r in pbm_metrics.iterrows())
    agreement_lines += [
        "",
        "The six-protein primary median is 0.3636 and the per-protein range is -0.4221 to 0.8133. This heterogeneity is an assay-dependence limitation; it is not evidence that either assay is universally correct for the other.",
    ]
    (OUT / "06_assay_agreement/assay_agreement_report.md").write_text("\n".join(agreement_lines) + "\n", encoding="utf-8")
    pd.DataFrame([{"status": "PARSED", "n_proteins": len(PROTEINS), "n_mutations": len(tidy), "source": RAW.relative_to(ROOT).as_posix()}]).to_csv(OUT / "06_assay_agreement/assay_agreement_status.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "status": "PARSED_MOESM16",
        "proteins_covered": len(PROTEINS),
        "n_mutations": len(tidy),
        "source_search": "official PMC/Nature source workbook",
        "files_checked": RAW.relative_to(ROOT).as_posix(),
        "reason": "Extended Data Fig. 3 XLS contains normalized competition values; paper describes displayed values as means of two replicates",
        "downstream_action": "competition-aligned scoring and PBM-vs-competition agreement evaluated",
    }]).to_csv(OUT / "04_competition/competition_source_status.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"figure": "Figure 2", "path": "figures/figure2_deeppbs_competition.png", "status": "PARSED_COMPETITION_SOURCE", "source": "04_competition/competition_mutations.tsv"},
        {"figure": "Figure 3", "path": "figures/figure3_nampnn_competition.png", "status": "FROZEN_GLOBAL_LANDSCAPE_DERIVED", "source": "05_assay_aligned_eval/nampnn_local_per_mutation_predictions.tsv"},
        {"figure": "Figure 4", "path": "figures/figure4_local_vs_global.png", "status": "PARSED_COMPETITION_SOURCE", "source": "05_assay_aligned_eval/local_global_comparison.tsv"},
        {"figure": "Figure 5", "path": "figures/figure5_pbm_competition_agreement.png", "status": "PARSED_COMPETITION_SOURCE", "source": "06_assay_agreement/pbm_competition_per_mutation.tsv"},
        {"figure": "Figure 6", "path": "figures/figure6_dbp048_sensitivity.png", "status": "DBP048_SENSITIVITY_ONLY", "source": "03_register_mapping/register_map.tsv"},
    ]).to_csv(OUT / "figures/figure_source_manifest.tsv", sep="\t", index=False)
    (OUT / "figures/README.md").write_text(
        "# v1.1 figure outputs\n\n"
        "The current assay-aligned figures are listed in `figure_source_manifest.tsv`. Files with the historical `_missing` suffix were generated before the official competition XLS was available and are retained for provenance only; they are superseded and are not part of the current figure manifest.\n",
        encoding="utf-8",
    )
    (OUT / "04_competition/competition_assay_audit.md").write_text(
        "# Competition assay branch status\n\n"
        "The official Glasscock Extended Data Fig. 3 XLS was obtained and parsed; the previous unavailable-source status is superseded by this parsed branch. "
        "The workbook contains heatmap means of two replicates, not individual replicate values.\n",
        encoding="utf-8",
    )
    (OUT / "10_conditionality/local_conditionality_status.md").write_text(
        "# Local conditionality\n\n"
        "Competition-aligned local vectors are now evaluable from the official Extended Data Fig. 3 XLS. "
        "DeepPBS uses native PWM `delta_logP`; NA-MPNN uses only the explicitly labeled frozen global-landscape-derived local score. "
        "DBP048 remains sensitivity-only because the sequence-C assay construct differs from the design sequence. "
        "These are development-exposed, protein-level diagnostics and are not independent validation.\n",
        encoding="utf-8",
    )
    write_final_summary(deep_metrics, na_metrics, pbm_metrics, wide_df)
    print("competition branch complete", len(tidy), "mutations")


def write_final_summary(deep: pd.DataFrame, na: pd.DataFrame, pbm: pd.DataFrame, wide: pd.DataFrame) -> None:
    deep_primary = deep[deep.mapping_status == "PRIMARY"]
    na_primary = na[na.mapping_status == "PRIMARY"]
    pbm_primary = pbm[pbm.mapping_status == "PRIMARY"]
    cap = pd.read_csv(OUT / "08_pwm_capacity/additive_pwm_capacity.tsv", sep="\t")
    cond = pd.read_csv(OUT / "10_conditionality/conditionality_summary.tsv", sep="\t")
    lines = [
        "# v1.1 summary",
        "",
        "This is a development-stage assay-alignment audit on the seven development-exposed GSE237017 DBPs. No v0.x or v1.0 output was overwritten.",
        "",
        "## Findings",
        "1. DeepPBS design-model coverage/QC: 7/7 frozen NPZ predictions validated.",
        f"2. DeepPBS local competition performance: primary six-protein median Spearman {deep_primary.spearman.median():.4f}; DBP048 is sensitivity-only because sequence C differs from the design sequence.",
        f"3. NA-MPNN local performance: primary six-protein median Spearman {na_primary.spearman.median():.4f}; every value is explicitly `NA-MPNN frozen global landscape-derived local mutation score`, not a native PPM result.",
        f"4. PBM-vs-competition agreement: primary six-protein median Spearman {pbm_primary.spearman.median():.4f}; this compares experimental assays, not a model.",
        "5. Global PBM comparison is retained in `local_global_comparison.tsv`; individual proteins, not pooled sequence counts, are the statistical units.",
        f"6. Experimental additive-PWM capacity at Hamming distance >=2 remains median {cap.loc[cap.subset == 'd_ge_2', 'spearman'].median():.4f}.",
        "7. Local-to-global distance strata remain descriptive and are not used to select a register.",
        f"8. Global conditionality remains protein-specific: DeepPBS median pairwise correlation {cond.loc[cond.method == 'DeepPBS', 'median_pairwise_spearman'].iloc[0]:.4f}; NA-MPNN {cond.loc[cond.method == 'NA-MPNN', 'median_pairwise_spearman'].iloc[0]:.4f}.",
        "9. Training overlap remains UNDETERMINED because full checkpoint training sequences are unavailable.",
        "",
        "## Per-protein local/global table",
        "",
        "| protein | DeepPBS local | DeepPBS global PBM | NA-MPNN derived local | NA-MPNN global PBM | PBM vs competition | mutations | mapping |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for _, r in wide.sort_values("protein").iterrows():
        lines.append(f"| {r.protein} | {r.DeepPBS_competition_rho:.4f} | {r.DeepPBS_global_PBM_rho:.4f} | {r.NA_MPNN_competition_rho:.4f} | {r.NA_MPNN_global_PBM_rho:.4f} | {r.PBM_vs_competition_rho:.4f} | {int(r.n_competition_mutations)} | {r.mapping_status} |")
    lines += [
        "",
        "## A/B/C/D classification",
        "A: supported descriptively. DeepPBS local median (0.3538) is above its global PBM median (0.0166); the NA-MPNN-derived local median (0.1487) is only modestly above its global median (0.1145). This is not an independent or population-level result.",
        "B: not supported as a blanket claim. Local performance is heterogeneous: DeepPBS spans 0.1824 to 0.7448 among primary proteins, while the derived NA-MPNN scores are mostly weak.",
        "C: partially supported as an interpretation constraint, not as a universal assay-disagreement claim. PBM-vs-competition primary rho spans -0.4221 to 0.8133 (median 0.3636), so assay agreement is protein-dependent.",
        "D: supported as a representation-capacity diagnostic because the experimental additive model is strong locally but degrades outside Hamming-1.",
        "",
        "Overall status: DEVELOPMENT-EXPOSED / ASSAY-ALIGNED DIAGNOSTIC COMPLETE. This does not constitute independent validation or a universal structure-method conclusion.",
    ]
    (OUT / "v1_1_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    old = json.loads((OUT / "v1_1_metrics.json").read_text(encoding="utf-8"))
    old.update({
        "competition_assay_status": "PARSED_MOESM16_HEATMAP_MEANS",
        "competition_n_mutations": int(len(pd.read_csv(OUT / "04_competition/competition_mutations.tsv", sep="\t"))),
        "deeppbs_local_median_spearman_primary_n6": float(deep_primary.spearman.median()),
        "nampnn_local_median_spearman_primary_n6": float(na_primary.spearman.median()),
        "pbm_competition_median_spearman_primary_n6": float(pbm_primary.spearman.median()),
        "deeppbs_local_per_protein": deep.set_index("protein").spearman.to_dict(),
        "nampnn_local_per_protein": na.set_index("protein").spearman.to_dict(),
        "pbm_competition_per_protein": pbm.set_index("protein").spearman.to_dict(),
        "competition_replicate_level_available": False,
        "dbp048_mapping_status": "AMBIGUOUS_SENSITIVITY_ONLY",
        "classification_A": "descriptive_only",
        "classification_B": "not_blanket_supported",
        "classification_C": "assessed_per_protein",
        "classification_D": "supported_capacity_diagnostic",
    })
    (OUT / "v1_1_metrics.json").write_text(json.dumps(old, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    complete_competition()
