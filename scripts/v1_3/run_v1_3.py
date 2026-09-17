from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
RAW_COMP = ROOT / "data/raw/v1_1_competition/source_data_extended_data_fig3.xls"
V11 = ROOT / "results/v1_1"
OUT = ROOT / "results/v1_3"
REPORT = ROOT / "reports/v1.3"
PROCESSED = ROOT / "data/processed"
OLD_PROTEINS = ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP035", "DBP048"]
PRIMARY_PROTEINS = ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP035"]
ALL_COMP = ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP023", "DBP035", "DBP048", "DBP056", "DBP062"]
LOCKED = ["DBP023", "DBP056", "DBP062"]
BASES = "ACGT"
SEED = 1301
N_BOOT = 2000
N_PERM = 2000


def safe_corr(x, y, kind="spearman"):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float((spearmanr if kind == "spearman" else pearsonr)(x, y).statistic)


def pairwise_accuracy(exp, pred):
    exp, pred = np.asarray(exp, float), np.asarray(pred, float)
    values = []
    for i in range(len(exp)):
        for j in range(i + 1, len(exp)):
            de, dp = exp[i] - exp[j], pred[i] - pred[j]
            if de == 0:
                continue
            values.append(1.0 if de * dp > 0 else 0.5 if dp == 0 else 0.0)
    return float(np.mean(values)) if values else np.nan


def ensure_dirs():
    for p in [OUT, REPORT, ROOT / "scripts/v1_3", PROCESSED]:
        p.mkdir(parents=True, exist_ok=True)


def parse_competition() -> pd.DataFrame:
    if not RAW_COMP.exists():
        raise FileNotFoundError(RAW_COMP)
    rows = []
    metadata = load_targets()
    for sheet in pd.ExcelFile(RAW_COMP).sheet_names:
        m = re.search(r"DBP(\d+)$", sheet)
        if not m:
            continue
        protein = f"DBP{int(m.group(1)):03d}"
        df = pd.read_excel(RAW_COMP, sheet_name=sheet)
        mut = df[df["position"].astype(str) != "WT"].copy()
        mut["position"] = mut["position"].astype(int)
        seq_by_pos = (mut.sort_values("position").drop_duplicates("position")
                      .set_index("position")["original_base"].astype(str).str.upper().to_dict())
        # The mutation workbook omits positions that were not assayed. For the
        # seven historical proteins, use the independently frozen assay target
        # metadata so an unmutated terminal base is not silently dropped.
        target = str(metadata.loc[protein, "experimental_assay_target"]) if protein in metadata.index else (structure_target_sequence(protein) or "".join(seq_by_pos[p] for p in sorted(seq_by_pos)))
        for i, r in mut.iterrows():
            wt, mb = str(r["original_base"]).upper(), str(r["new_base"]).upper()
            raw = float(r["Median PE/FITC (Normalized)"])
            rows.append({
                "protein_id": protein, "assay": "competition", "target_sequence": target,
                "position": int(r["position"]), "wt_base": wt, "mutant_base": mb,
                "replicate": "mean_of_two_replicates", "raw_measurement": raw,
                "normalized_effect": -raw,
                "source_file": RAW_COMP.relative_to(ROOT).as_posix(),
                "source_sheet": sheet, "source_cell": f"A{i + 2}:E{i + 2}",
                "development_status": "development_exposed" if protein in OLD_PROTEINS else "locked_holdout",
                "holdout_status": "DEVELOPMENT_EXPOSED" if protein in OLD_PROTEINS else "LOCKED_HOLDOUT",
                "replicate_level_available": False,
            })
    out = pd.DataFrame(rows).sort_values(["protein_id", "position", "mutant_base"]).reset_index(drop=True)
    expected = set(ALL_COMP)
    assert set(out.protein_id) == expected, (set(out.protein_id), expected)
    assert not out.duplicated(["protein_id", "position", "wt_base", "mutant_base"]).any()
    for (protein, pos), g in out.groupby(["protein_id", "position"]):
        assert len(g) == 3 and len(set(g.mutant_base)) == 3 and all(b in BASES for b in g.mutant_base)
    out.to_parquet(PROCESSED / "v1_3_competition_mutations.parquet", index=False)
    # A consensus file is intentionally identical to the published two-replicate mean;
    # the raw table above remains the canonical tidy table and is never averaged again.
    out.to_parquet(PROCESSED / "v1_3_competition_mutations_consensus.parquet", index=False)
    return out


def load_targets():
    p = ROOT / "metadata/v0_3_1/designed_dbp_target_definitions.csv"
    d = pd.read_csv(p)
    d["protein_id"] = d.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}")
    return d.set_index("protein_id")


def structure_status(protein):
    p = ROOT / "data/raw/v1_0_structure/design_pdbs/design_pdbs" / f"{protein}.pdb"
    return ("available_file" if p.exists() else "not_found_in_checked_out_structure_dir", p.relative_to(ROOT).as_posix() if p.exists() else "")


def structure_target_sequence(protein):
    """Extract the first strand from a bundled DNA duplex when available."""
    p = ROOT / "data/raw/v1_0_structure/design_pdbs/design_pdbs" / f"{protein}.pdb"
    if not p.exists():
        return ""
    code = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "A": "A", "C": "C", "G": "G", "T": "T"}
    seen, keys = [], set()
    for line in p.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        residue, chain, number = line[17:20].strip(), line[21].strip(), line[22:26].strip()
        key = (chain, number)
        if chain == "B" and residue in code and key not in keys:
            keys.add(key)
            seen.append(code[residue])
    duplex = "".join(seen)
    return duplex[: len(duplex) // 2] if len(duplex) >= 2 and len(duplex) % 2 == 0 else duplex


def write_audit(comp):
    targets = load_targets()
    pbm = pd.read_parquet(ROOT / "data/processed/v0_3_1/designed_dbp_upbm_rc_class_v0_3_1.parquet")
    pbm["protein_id"] = pbm.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}")
    deep = pd.read_csv(V11 / "05_assay_aligned_eval/deeppbs_per_mutation_predictions.tsv", sep="\t")
    rows = []
    for protein in ALL_COMP:
        g = comp[comp.protein_id == protein]
        positions, variants = g.position.nunique(), len(g)
        structure, structure_file = structure_status(protein)
        has_pbm = protein in set(pbm.protein_id)
        old = protein in OLD_PROTEINS
        analysis_status = "LOCKED_HOLDOUT" if protein in LOCKED else ("PRIMARY" if protein in PRIMARY_PROTEINS else "SENSITIVITY_ONLY")
        target = str(g.target_sequence.iloc[0])
        target_source = "metadata/v0_3_1 experimental_assay_target" if protein in targets.index else "bundled PDB DNA chain B first strand"
        rows.append({
            "protein_id": protein, "assay": "competition", "target_sequence": target,
            "mutation_coverage": "3 substitutions per position", "n_positions": positions,
            "n_variants": variants, "replicates": "two-replicate mean only; raw replicates unavailable",
            "uPBM": "available (2 replicates)" if has_pbm else "not available in local GSE237017 benchmark",
            "competition": "available; XLS mean of two replicates",
            "structure_availability": structure, "structure_file": structure_file,
            "development_status": "development_exposed" if old else "not used in v0.x-v1.2 outputs found by audit",
            "holdout_status": "DEVELOPMENT_EXPOSED" if old else "LOCKED_HOLDOUT",
            "analysis_status": analysis_status,
            "deepPBS_prediction": "available" if protein in set(deep.protein) else "not available",
            "NA_MPNN_prediction": "available local landscape only" if protein in set(deep.protein) else "not available",
            "data_source": RAW_COMP.relative_to(ROOT).as_posix(), "target_sequence_source": target_source,
            "missing_items": "individual competition replicates; DBP35opt pair" + ("; PBM" if not has_pbm else ""),
        })
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "data_audit.tsv", sep="\t", index=False)
    lines = [
        "# v1.3 data audit", "",
        "Audit date: 2026-09-17. Counts below were computed from the checked-out files; no quantities were inferred.", "",
        "## Scope and status", "",
        "The competition source is `data/raw/v1_1_competition/source_data_extended_data_fig3.xls` (10 sheets). The seven proteins used throughout v0.x-v1.2 are development-exposed. DBP023, DBP056, and DBP062 occur in the source workbook and in the raw design-PDB directory but were absent from prior v1.1/v1.2 prediction and metric tables found by this audit; they are frozen as `LOCKED_HOLDOUT` before v1.3 method results. The primary method summary uses DBP001/003/005/006/009/035; DBP048 remains `SENSITIVITY_ONLY` because its mapping ambiguity was documented before v1.3.", "",
        "The workbook contains a normalized PE/FITC value described as the mean of two replicates. It does not contain the two underlying replicate values. Therefore the v1.3 tidy table preserves one row per published mean and does not manufacture replicate rows. PBM has two replicate-level source files for the seven GSE237017 proteins.", "",
        "`DBP35opt` was not found in repository paths, metadata, cached predictions, or reports searched for the exact identifiers `DBP35opt`/`DBP035opt`; the DBP35 -> DBP35opt case study is unavailable until a source is supplied.", "",
        "## Per-protein inventory", "",
        "| protein | assay | target sequence | target source | positions | variants | replicates | structure | analysis status | development/holdout | DeepPBS | missing |",
        "|---|---|---|---|---:|---:|---|---|---|---|---|---|",
    ]
    for _, r in audit.iterrows():
        lines.append(f"| {r.protein_id} | {r.assay} | `{r.target_sequence}` | {r.target_sequence_source} | {r.n_positions} | {r.n_variants} | {r.replicates} | {r.structure_availability} | {r.analysis_status} | {r.holdout_status} | {r.deepPBS_prediction} | {r.missing_items} |")
    lines += [
        "", "## Existing development artifacts", "",
        "- uPBM: `data/processed/v0_3_1/designed_dbp_upbm_rc_class_v0_3_1.parquet`, 57,344 protein-RC-class rows, seven proteins, two PBM replicates per protein.",
        "- Competition: `data/raw/v1_1_competition/source_data_extended_data_fig3.xls`, parsed into `data/processed/v1_3_competition_mutations.parquet`.",
        "- DeepPBS: frozen NPZ bundle and v1.1 per-mutation table for seven proteins; DBP048 is mapping sensitivity-only in prior reports.",
        "- Conditionality: frozen v0.6 M3 protein- and target-shuffle tables exist under `results/v0_6_representation/`; v1.3 standardizes them without rerunning predictions.",
        "- Existing position/local metrics: `results/v1_2/04_position_decomposition/`, `results/v1_2/07_assay_triangle/`, and `results/v1_2/09_local_conditionality/`; v1.3 recomputes metrics from the tidy table and frozen predictions without changing their definitions after contract freeze.",
        "- No DBP023/056/062 PBM or model prediction was found; these proteins are not silently promoted to an accuracy benchmark.",
        "",
        "## Audit limitations", "",
        "Structure-file presence is reported separately from validated structure/assay alignment. The three locked holdouts have raw PDB files but no checked-out v0.x-v1.2 structure metadata rows. Replicate reproducibility for competition cannot be estimated from this workbook; only PBM replicate data are available locally.",
    ]
    (REPORT / "V1_3_DATA_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return audit


def write_contract():
    text = """# V1.3 Evaluation Contract

Status: **FROZEN 2026-09-17**. This contract is written before v1.3 results. Metric definitions may change only for a demonstrated implementation bug; a result that is inconvenient is not a bug.

## Analysis units and independence

The primary statistical unit is the protein. Mutation rows are repeated observations within a protein landscape, not independent biological replicates. All primary summaries report per-protein values, the number of proteins, median, mean as an auxiliary summary, protein-level bootstrap 95% intervals, and leave-one-protein-out values where implemented. Development-exposed proteins are not independent validation.

The seven GSE237017 proteins and all v1.1/v1.2 model outputs are `development_exposed`. DBP023, DBP056, and DBP062 are frozen `LOCKED_HOLDOUT` because the source workbook contains them but the audit found no prior prediction/metric artifact using them. They cannot be used for tuning. This status is frozen before reading v1.3 benchmark results.

The primary method summary contains DBP001, DBP003, DBP005, DBP006, DBP009, and DBP035. DBP048 remains `SENSITIVITY_ONLY` because v1.1/v1.2 had already documented an ambiguous local mapping; this exclusion was inherited before v1.3 outcomes and is not result-driven. An all-seven descriptive summary is retained and explicitly labeled.

## Measurement direction and ties

Competition raw measurement is normalized PE/FITC. The frozen monotone effect is `normalized_effect = -raw_measurement`; larger values mean stronger competition. The raw value remains in the tidy table. PBM and competition are different assay estimands and are never treated as interchangeable ground truth. Spearman uses average ranks (the scipy default). Pairwise identity accuracy gives a comparable pair score of 1 for a correct strict order, 0 for an incorrect strict order, and 0.5 for a prediction tie; experimental ties are excluded from the denominator. Its chance baseline is 0.5.

## A. Global mutation-effect endpoint

For every protein and frozen method, calculate Spearman and Pearson over all available mutation rows for that protein. Report the per-protein values. A descriptive rank-based 95% interval is obtained by deterministic position-stratified bootstrap of mutation rows within that protein; it is not a biological-replicate interval. The primary across-protein interval is a protein bootstrap (2,000 resamples, seed 1301) of the per-protein statistic. No pooled mutation-row p-value is primary.

## B. Position sensitivity endpoint

For protein p and position j, the primary position sensitivity is `mean_abs_effect(j) = mean_b(abs(effect(j,b)))` over the three substitutions. This definition is fixed because competition direction is an assay-specific signed signal and absolute magnitude directly represents sensitivity without choosing a deleterious sign threshold. A WT-relative "mean deleterious effect" is not identifiable from this workbook because it provides normalized competitor signal rather than a WT-relative signed change. `mean_signed_effect` and `max_abs_effect` are retained as secondary diagnostics. Predicted values use the same transformation. The endpoint is Spearman between experimental and predicted position summaries, per protein, with protein-level bootstrap intervals.

## C. Within-position nucleotide identity endpoint

For each position, residualize both vectors by their within-position mean: `delta(j,b) = effect(j,b) - mean_b(effect(j,b))`. Primary metrics are pooled within-position residual Spearman and Pearson within each protein, plus mean pairwise accuracy across positions. Identity pair comparisons are only within a position. Position-centering removes position sensitivity; no global centering is used. Chance is 0 for residual correlation under a symmetric random-order null and 0.5 for pairwise accuracy. The reported Kendall-style concordance is the equivalent transform `2 * pairwise_accuracy - 1`, with chance 0; it is not treated as an additional independent endpoint. Within-position base-label permutation shuffles predicted mutant labels within each position (2,000 permutations, seed 1301); this is the primary identity null.

## D. Protein conditionality

For cached global landscapes, report pairwise protein Spearman, mean absolute prediction shift between proteins, and per-7-mer variance across protein identities. The landscape conditionality effect size is the mean across-7-mer between-protein variance divided by total variance across all protein/7-mer scores; it is descriptive and is not an accuracy claim. For the existing frozen v0.6 M3 protein- and target-shuffle predictions, preserve prediction correlation, mean absolute score shift, and the historical normalized rank-shift definition `1 - prediction_correlation`. These are development-exposed diagnostics, not a new holdout evaluation. Protein-label permutation is at the protein unit. If a future method has no target-conditioned artifact, target shuffle is reported unavailable rather than simulated.

The real perturbation endpoint requires matched DBP35 and DBP35opt experimental and predicted landscapes. The audit found no DBP35opt source, so this endpoint is pre-specified but currently `NOT EVALUABLE`, not replaced by a synthetic proxy.

## Nulls and uncertainty

Position-label permutation shuffles predicted position summaries across positions within protein. Within-position base-label permutation shuffles predicted mutant identity within each position. Protein-label permutation is applied only at the protein unit. Null distributions and observed values are written to source tables. Correlation endpoints use two-sided empirical permutation p-values around their chance value of zero; pairwise accuracy uses the upper tail around chance 0.5. The add-one correction is used. Leave-one-protein-out sensitivity is descriptive; no protein is selected after seeing its result.

## Method inclusion

The first benchmark includes frozen DeepPBS, frozen NA-MPNN landscape-derived scores when the local register has coverage, experimental PBM as an assay comparator, and position-only summaries. No new neural model is trained in v1.3. Methods lacking a complete, auditable prediction are marked unavailable in the method reproducibility record and are not backfilled.

## Figure 2 data contract

Figure 2 source data must contain one row per method/protein with global Spearman, position-sensitivity Spearman, within-position residual Spearman, within-position pairwise accuracy, mutation count, holdout/development status, and the fixed chance baselines. The figure is generated only from this source table; no hand-edited numbers are permitted.
"""
    (ROOT / "V1_3_EVALUATION_CONTRACT.md").write_text(text, encoding="utf-8")


def position_metrics(g, pred_col):
    exp_pos, pred_pos = [], []
    exp_signed, pred_signed, exp_max_abs, pred_max_abs = [], [], [], []
    residual_exp, residual_pred, pair = [], [], []
    for pos, q in g.groupby("position", sort=True):
        e, p = q.normalized_effect.to_numpy(float), q[pred_col].to_numpy(float)
        exp_pos.append(np.mean(np.abs(e))); pred_pos.append(np.mean(np.abs(p)))
        exp_signed.append(np.mean(e)); pred_signed.append(np.mean(p))
        exp_max_abs.append(np.max(np.abs(e))); pred_max_abs.append(np.max(np.abs(p)))
        residual_exp.extend(e - np.mean(e)); residual_pred.extend(p - np.mean(p))
        pair.append(pairwise_accuracy(e - np.mean(e), p - np.mean(p)))
    pairwise = float(np.nanmean(pair))
    return {
        "position_spearman": safe_corr(exp_pos, pred_pos),
        "position_pearson": safe_corr(exp_pos, pred_pos, "pearson"),
        "position_mean_signed_spearman": safe_corr(exp_signed, pred_signed),
        "position_mean_signed_pearson": safe_corr(exp_signed, pred_signed, "pearson"),
        "position_max_abs_spearman": safe_corr(exp_max_abs, pred_max_abs),
        "position_max_abs_pearson": safe_corr(exp_max_abs, pred_max_abs, "pearson"),
        "within_residual_spearman": safe_corr(residual_exp, residual_pred),
        "within_residual_pearson": safe_corr(residual_exp, residual_pred, "pearson"),
        "within_pairwise_accuracy": pairwise,
        "within_kendall_style_concordance": 2 * pairwise - 1,
        "n_positions": len(exp_pos),
    }


def prediction_sources():
    """Load the frozen local prediction tables without running model inference."""
    deep = pd.read_csv(V11 / "05_assay_aligned_eval/deeppbs_per_mutation_predictions.tsv", sep="\t")
    deep = deep.rename(columns={"predicted_delta_logP": "predicted"})
    pbm = pd.read_csv(V11 / "06_assay_agreement/pbm_competition_per_mutation.tsv", sep="\t")
    pbm = pbm.rename(columns={"predicted_delta": "predicted"})
    na = pd.read_csv(V11 / "05_assay_aligned_eval/nampnn_local_per_mutation_predictions.tsv", sep="\t")
    na = na.rename(columns={"predicted_delta": "predicted"})
    return {"DeepPBS": deep, "PBM_experimental": pbm, "NA-MPNN": na}


def join_competition_predictions(pred, protein, comp):
    q = pred[pred.protein == protein].copy()
    if q.empty:
        return q
    q = q.merge(
        comp[["protein_id", "position", "wt_base", "mutant_base", "normalized_effect", "holdout_status"]],
        left_on=["protein", "assay_position", "wt_base", "mut_base"],
        right_on=["protein_id", "position", "wt_base", "mutant_base"],
        how="inner",
    )
    return q.reset_index(drop=True)


def stratified_rank_ci(g, pred_col, kind="spearman", rng=None):
    rng = np.random.default_rng(0) if rng is None else rng
    vals = []
    groups = [q for _, q in g.groupby("position")]
    for _ in range(500):
        rows = []
        for q in groups:
            rows.append(q.iloc[rng.integers(0, len(q), len(q))])
        b = pd.concat(rows, ignore_index=True)
        vals.append(safe_corr(b.normalized_effect, b[pred_col], kind))
    vals = np.asarray(vals, float); vals = vals[np.isfinite(vals)]
    return (float(np.quantile(vals, .025)), float(np.quantile(vals, .975))) if len(vals) else (np.nan, np.nan)


def per_protein_decomposition(comp):
    rows = []
    for method, pred in prediction_sources().items():
        for protein in OLD_PROTEINS:
            q = join_competition_predictions(pred, protein, comp)
            if q.empty:
                continue
            global_s = safe_corr(q.normalized_effect, q.predicted)
            global_p = safe_corr(q.normalized_effect, q.predicted, "pearson")
            ci_lo, ci_hi = stratified_rank_ci(q, "predicted", rng=np.random.default_rng(SEED + OLD_PROTEINS.index(protein)))
            m = position_metrics(q, "predicted")
            mapping = ";".join(sorted(q.mapping_status.dropna().astype(str).unique())) if "mapping_status" in q else "not_recorded"
            rows.append({
                "method": method,
                "protein_id": protein,
                "analysis_status": "PRIMARY" if protein in PRIMARY_PROTEINS else "SENSITIVITY_ONLY",
                "mapping_status": mapping,
                "coverage_scope": "full_competition_landscape" if method == "DeepPBS" else "fixed_7mer_window",
                "n_mutations": len(q),
                "holdout_status": "DEVELOPMENT_EXPOSED",
                "global_spearman": global_s,
                "global_pearson": global_p,
                "global_rank_ci_low": ci_lo,
                "global_rank_ci_high": ci_hi,
                **m,
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "decomposition_per_protein.tsv", sep="\t", index=False)
    return out


def protein_bootstrap_summary(per):
    metrics = ["global_spearman", "position_spearman", "within_residual_spearman", "within_pairwise_accuracy"]
    rows = []
    rng = np.random.default_rng(SEED)
    analysis_sets = {
        "primary_six": per[per.analysis_status == "PRIMARY"],
        "all_seven_descriptive": per,
    }
    for analysis_set, table in analysis_sets.items():
        for method, g in table.groupby("method"):
            for metric in metrics:
                x = g[metric].dropna().to_numpy(float)
                draws = np.array([np.median(rng.choice(x, len(x), replace=True)) for _ in range(N_BOOT)]) if len(x) else np.array([])
                rows.append({"analysis_set": analysis_set, "method": method, "metric": metric, "n_proteins": len(x), "median": np.median(x) if len(x) else np.nan,
                             "mean": np.mean(x) if len(x) else np.nan,
                             "ci_low": np.quantile(draws, .025) if len(draws) else np.nan,
                             "ci_high": np.quantile(draws, .975) if len(draws) else np.nan,
                             "chance_baseline": .5 if metric == "within_pairwise_accuracy" else 0.0})
    out = pd.DataFrame(rows); out.to_csv(OUT / "decomposition_summary.tsv", sep="\t", index=False); return out


def leave_one_protein_out(per):
    metrics = ["global_spearman", "position_spearman", "within_residual_spearman", "within_pairwise_accuracy"]
    rows = []
    primary = per[per.analysis_status == "PRIMARY"]
    for method, g in primary.groupby("method"):
        for omitted in sorted(g.protein_id.unique()):
            remaining = g[g.protein_id != omitted]
            for metric in metrics:
                x = remaining[metric].dropna().to_numpy(float)
                rows.append({
                    "analysis_set": "primary_six",
                    "method": method,
                    "omitted_protein": omitted,
                    "metric": metric,
                    "n_proteins": len(x),
                    "median": float(np.median(x)) if len(x) else np.nan,
                    "mean": float(np.mean(x)) if len(x) else np.nan,
                })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "leave_one_protein_out.tsv", sep="\t", index=False)
    return out


def permutation_nulls(comp, per):
    summary_rows, distribution_rows = [], []
    for method_index, (method, pred) in enumerate(prediction_sources().items()):
        for protein_index, protein in enumerate(OLD_PROTEINS):
            q = join_competition_predictions(pred, protein, comp)
            if q.empty:
                continue
            obs = position_metrics(q, "predicted")
            rng = np.random.default_rng(SEED + 100 * method_index + protein_index)
            exp = q.normalized_effect.to_numpy(float)
            predicted = q.predicted.to_numpy(float)
            groups = [x.index.to_numpy() for _, x in q.groupby("position", sort=True)]
            exp_pos = np.array([np.mean(np.abs(exp[idx])) for idx in groups])
            pred_pos = np.array([np.mean(np.abs(predicted[idx])) for idx in groups])
            exp_residual = np.concatenate([exp[idx] - np.mean(exp[idx]) for idx in groups])
            null_values = {
                "position_spearman": [],
                "within_residual_spearman": [],
                "within_pairwise_accuracy": [],
            }
            for permutation in range(N_PERM):
                null_values["position_spearman"].append(safe_corr(exp_pos, rng.permutation(pred_pos)))
                shuffled = predicted.copy()
                for idx in groups:
                    shuffled[idx] = rng.permutation(shuffled[idx])
                pred_residual = np.concatenate([shuffled[idx] - np.mean(shuffled[idx]) for idx in groups])
                null_values["within_residual_spearman"].append(safe_corr(exp_residual, pred_residual))
                null_values["within_pairwise_accuracy"].append(
                    float(np.nanmean([pairwise_accuracy(exp[idx], shuffled[idx]) for idx in groups]))
                )
            definitions = [
                ("position_label", "position_spearman", obs["position_spearman"], "two_sided", 0.0),
                ("within_position_base_label", "within_residual_spearman", obs["within_residual_spearman"], "two_sided", 0.0),
                ("within_position_base_label", "within_pairwise_accuracy", obs["within_pairwise_accuracy"], "upper", 0.5),
            ]
            for null_type, metric, observed, tail, chance in definitions:
                values = np.asarray(null_values[metric], float)
                values = values[np.isfinite(values)]
                if tail == "two_sided":
                    count = np.sum(np.abs(values - chance) >= abs(observed - chance))
                else:
                    count = np.sum(values >= observed)
                p_value = float((1 + count) / (1 + len(values)))
                analysis_status = "PRIMARY" if protein in PRIMARY_PROTEINS else "SENSITIVITY_ONLY"
                summary_rows.append({
                    "method": method,
                    "protein_id": protein,
                    "analysis_status": analysis_status,
                    "null_type": null_type,
                    "metric": metric,
                    "observed": observed,
                    "chance_baseline": chance,
                    "null_median": np.nanmedian(values),
                    "null_q025": np.nanquantile(values, .025),
                    "null_q975": np.nanquantile(values, .975),
                    "empirical_p": p_value,
                    "tail": tail,
                    "n_permutations": len(values),
                })
                distribution_rows.extend({
                    "method": method,
                    "protein_id": protein,
                    "analysis_status": analysis_status,
                    "null_type": null_type,
                    "metric": metric,
                    "permutation": i,
                    "null_value": value,
                } for i, value in enumerate(values))
    out = pd.DataFrame(summary_rows)
    out.to_csv(OUT / "permutation_nulls.tsv", sep="\t", index=False)
    pd.DataFrame(distribution_rows).to_csv(OUT / "permutation_distributions.tsv", sep="\t", index=False)
    return out


def conditionality():
    p = V11 / "09_global_pwm_projection/unified_global_landscapes.tsv"
    if not p.exists(): return pd.DataFrame()
    d = pd.read_csv(p, sep="\t")
    rows = []
    for method, g in d.groupby("method"):
        wide = g.pivot(index="canonical_7mer", columns="protein", values="predicted_score")
        cols = [x for x in OLD_PROTEINS if x in wide.columns]
        arr = wide[cols].to_numpy(float)
        between = np.nanvar(arr, axis=1)
        total = np.nanvar(arr)
        for a in range(len(cols)):
            for b in range(a + 1, len(cols)):
                rows.append({"method": method, "protein_a": cols[a], "protein_b": cols[b], "spearman": safe_corr(arr[:, a], arr[:, b]), "mean_abs_prediction_shift": float(np.nanmean(np.abs(arr[:, a] - arr[:, b]))), "n_7mers": len(wide)})
        rows.append({"method": method, "protein_a": "SUMMARY", "protein_b": "SUMMARY", "spearman": np.nan, "mean_abs_prediction_shift": float(np.nanmean(between)), "n_7mers": len(wide), "protein_variance_fraction": float(np.nanmean(between) / total) if total else np.nan})
    out = pd.DataFrame(rows); out.to_csv(OUT / "conditionality_global.tsv", sep="\t", index=False); return out


def cached_shuffle_diagnostics():
    """Standardize existing v0.6 shuffle caches without rerunning predictions."""
    source_dir = ROOT / "results/v0_6_representation"
    rows = []
    definitions = [
        ("protein", "protein_shuffle.csv", "replacement_protein_id", "protein_conditioning_effect_size"),
        ("target", "target_shuffle.csv", "replacement_target_dbp_id", "target_conditioning_effect_size"),
    ]
    for shuffle_type, filename, replacement_col, effect_col in definitions:
        table = pd.read_csv(source_dir / filename)
        for _, row in table.iterrows():
            rows.append({
                "model": "M3",
                "shuffle_type": shuffle_type,
                "split_type": row.split_type,
                "fold_id": row.fold_id,
                "seed": int(row.seed),
                "protein_id": f"DBP{int(str(row.dbp_id).replace('DBP', '')):03d}",
                "replacement_id": f"DBP{int(str(row[replacement_col]).replace('DBP', '')):03d}",
                "prediction_correlation": row.prediction_correlation,
                "mean_abs_prediction_shift": row.mean_abs_score_change,
                "normalized_rank_shift": row[effect_col],
                "retrained": bool(row.retrained),
                "development_status": "development_exposed_historical_v0_6",
                "source_file": (source_dir / filename).relative_to(ROOT).as_posix(),
            })
    detail = pd.DataFrame(rows)
    detail.to_csv(OUT / "cached_shuffle_diagnostics.tsv", sep="\t", index=False)
    summary = pd.read_csv(source_dir / "conditionality_diagnostics.csv")
    summary["source_file"] = "results/v0_6_representation/conditionality_diagnostics.csv"
    summary["development_status"] = "development_exposed_historical_v0_6"
    summary.to_csv(OUT / "cached_shuffle_summary.tsv", sep="\t", index=False)
    return detail, summary


def figure2(per):
    rows = []
    for _, r in per.iterrows():
        rows.append({"method": r.method, "protein_id": r.protein_id, "global_spearman": r.global_spearman, "position_spearman": r.position_spearman,
                     "within_position_residual_spearman": r.within_residual_spearman, "within_position_pairwise_accuracy": r.within_pairwise_accuracy,
                     "within_kendall_style_concordance": r.within_kendall_style_concordance,
                     "n_mutations": r.n_mutations, "n_positions": r.n_positions, "holdout_status": r.holdout_status,
                     "analysis_status": r.analysis_status, "mapping_status": r.mapping_status,
                     "coverage_scope": r.coverage_scope, "global_chance": 0.0, "position_chance": 0.0,
                     "within_residual_chance": 0.0, "within_pairwise_chance": .5,
                     "within_kendall_chance": 0.0})
    out = pd.DataFrame(rows); out.to_csv(OUT / "figure2_data.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    labels = [("global_spearman", "Global mutation effect"), ("position_spearman", "Position sensitivity"), ("within_position_residual_spearman", "Within-position identity")]
    methods = ["DeepPBS", "NA-MPNN", "PBM_experimental"]
    colors = dict(zip(methods, plt.rcParams["axes.prop_cycle"].by_key()["color"][:3]))
    for ax, (col, title) in zip(axes, labels):
        for method_index, method in enumerate(methods):
            g = out[out.method == method]
            primary = g[g.analysis_status == "PRIMARY"].sort_values("protein_id")
            sensitivity = g[g.analysis_status == "SENSITIVITY_ONLY"]
            offsets = np.linspace(-0.11, 0.11, len(primary))
            ax.scatter(method_index + offsets, primary[col], color=colors[method], s=38)
            ax.scatter([method_index] * len(sensitivity), sensitivity[col], marker="x", color="black", s=48, linewidths=1.2)
        ax.axhline(0, color="0.75", lw=.8); ax.set_title(title); ax.tick_params(axis="x", rotation=45)
        ax.set_xticks(range(len(methods)), methods)
    axes[0].set_ylabel("per-protein Spearman")
    handles = [Line2D([0], [0], marker="o", linestyle="none", color=colors[m], label=m, markersize=7) for m in methods]
    handles.append(Line2D([0], [0], marker="x", linestyle="none", color="black", label="DBP048 sensitivity-only", markersize=7))
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, .86)); fig.savefig(OUT / "figure2_global_position_identity.png", dpi=220); plt.close(fig)


def update_progress(per, audit):
    summary = pd.read_csv(OUT / "decomposition_summary.tsv", sep="\t")
    lines = ["# v1.3 progress", "", "## Completed", "", "- Audited README, PROGRESS, v1.1/v1.2 reports, scripts, results, metadata, raw competition workbook, PBM, structure files, and cached DeepPBS/NA-MPNN artifacts.", "- Frozen `V1_3_DATA_AUDIT.md` and root `V1_3_EVALUATION_CONTRACT.md` before v1.3 benchmark results.", "- Parsed all 10 competition sheets into `data/processed/v1_3_competition_mutations.parquet` and a consensus copy. Published values are two-replicate means; individual replicate rows are unavailable.", "- Recomputed global, position-sensitivity, and within-position residual identity metrics for frozen DeepPBS/PBM/NA-MPNN inputs. Generated deterministic permutation null distributions and empirical p-values, protein bootstrap summaries, leave-one-protein-out sensitivity, conditionality diagnostics, and Figure 2 source data/PNG.", "- Standardized the existing v0.6 M3 protein- and target-shuffle caches without rerunning inference; these remain development-exposed historical diagnostics.", "", "## Key audit numbers", "", f"- Competition: {len(audit)} proteins, {int(audit.n_variants.sum())} mutation rows, {int(audit.n_positions.sum())} position records.", f"- Development-exposed: {sum(audit.holdout_status == 'DEVELOPMENT_EXPOSED')}; locked holdout: {sum(audit.holdout_status == 'LOCKED_HOLDOUT')} (DBP023/056/062).", "- Primary reported method set: DBP001/003/005/006/009/035. DBP048 remains sensitivity-only because its local mapping ambiguity predates v1.3.", f"- DBP35opt: not found in checked-out data/metadata/results; real perturbation endpoint remains NOT EVALUABLE."]
    primary_summary = summary[summary.analysis_set == "primary_six"]
    for method, g in primary_summary.groupby("method"):
        x = g.set_index("metric"); lines.append(f"- {method}: global median {x.loc['global_spearman','median']:.3f}; position median {x.loc['position_spearman','median']:.3f}; within-identity residual median {x.loc['within_residual_spearman','median']:.3f}; pairwise median {x.loc['within_pairwise_accuracy','median']:.3f} (chance 0.5).")
    lines += ["", "## Stable conclusions", "", "- The benchmark now distinguishes global ranking, position sensitivity, and within-position identity without pooling mutation rows as independent proteins.", "- In the primary six, DeepPBS position-sensitivity agreement (median rho 0.406) is stronger than residual identity agreement (median rho 0.114); pairwise identity accuracy is 0.506, effectively the 0.5 chance baseline at the protein-summary level.", "- The seven historical proteins remain development-exposed; DBP023/056/062 are not yet method-evaluable despite having competition rows.", "- Cross-assay and competition-replicate variance decomposition is limited by absent individual competition replicate values; PBM replicate information is available.", "", "## Not yet concluded", "", "- DeepPBS DBP035 is the sole primary protein with nominal within-position evidence in both residual rho (0.444, permutation p=0.027) and pairwise accuracy (0.667, p=0.029). It is development-exposed and these p-values are not multiplicity-adjusted, so this is a case-level lead rather than validation.", "- No claim is made about DBP35 -> DBP35opt because the paired protein is absent.", "- No new GNN/model development is justified before reviewing the frozen decomposition and holdout availability.", "", "## Next single priority", "", "Obtain an auditable frozen prediction for DBP023/056/062 without using their competition values for model or metric tuning; otherwise the benchmark cannot make a prospective holdout claim."]
    (REPORT / "V1_3_PROGRESS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    old = ROOT / "PROGRESS.md"
    marker = "\n\n### v1.3 Conditionality & Local Effect Benchmark (2026-09-17)"
    existing = old.read_text(encoding="utf-8")
    if marker in existing:
        existing = existing.split(marker, 1)[0].rstrip()
    block = "\n\n### v1.3 Conditionality & Local Effect Benchmark (2026-09-17)\n\n- Audited the repository and froze `reports/v1.3/V1_3_DATA_AUDIT.md` plus `V1_3_EVALUATION_CONTRACT.md`.\n- Parsed the 10-sheet competition workbook into the v1.3 tidy/consensus Parquet tables; DBP001/003/005/006/009/035/048 are development-exposed and DBP023/056/062 are locked holdouts.\n- Recomputed global, position-sensitivity, and within-position identity metrics from frozen v1.1 inputs; generated full permutation nulls, empirical p-values, protein-bootstrap and leave-one-protein-out summaries, and `results/v1_3/figure2_data.tsv`. The primary set is six proteins; DBP048 remains sensitivity-only.\n- DBP35opt and individual competition replicate values remain unavailable; no new model was trained. See `reports/v1.3/V1_3_PROGRESS.md`.\n"
    old.write_text(existing + block, encoding="utf-8")


def method_reproducibility():
    text = """# v1.3 method reproducibility

This record freezes method inclusion before interpreting Figure 2. A method is included in the local decomposition only when a cached prediction can be joined to the competition key `(protein, position, wt_base, mutant_base)` without refitting or selecting a register from the v1.3 outcome.

| method | local input | status | reason |
|---|---|---|---|
| DeepPBS | `results/v1_1/05_assay_aligned_eval/deeppbs_per_mutation_predictions.tsv` | INCLUDED | Frozen native PWM delta-logP; all seven development-exposed proteins, fixed v1.1 register. |
| NA-MPNN | `results/v1_1/05_assay_aligned_eval/nampnn_local_per_mutation_predictions.tsv` | INCLUDED_WITH_COVERAGE_LIMIT | Frozen global-landscape-derived local scores; only the 7 positions covered by the fixed mapped window (21 mutations/protein). No new inference. |
| PBM experimental comparator | `results/v1_1/06_assay_agreement/pbm_competition_per_mutation.tsv` | INCLUDED_AS_ASSAY_COMPARATOR | Experimental PBM delta from the frozen 7-mer register; not a computational method and not interchangeable ground truth. |
| M0/M1/M1c/M2/M3 | historical v0.5-v0.6 global 7-mer artifacts | CONDITIONALITY_CACHE_ONLY | Available artifacts do not expose a frozen competition-key prediction table for Figure 2. Existing M3 protein/target shuffle tables are standardized in `results/v1_3/cached_shuffle_diagnostics.tsv` without new inference. |
| historical sequence baseline/SimplePC | v0.3-v0.4 benchmark artifacts | NOT_INCLUDED_IN_LOCAL_FIGURE_2 | Sequence/PBM benchmark outputs are not a frozen per-mutation prediction table for this competition assay and cannot be silently treated as a local model. |
| DBP023/DBP056/DBP062 | competition source only | LOCKED_HOLDOUT_UNEVALUATED | No cached PBM, DeepPBS, NA-MPNN, M0-M3, or sequence baseline prediction was found for these proteins. They remain untouched for future prospective evaluation. |
| DBP35 -> DBP35opt | no paired source | NOT_EVALUABLE | Exact DBP35opt identifiers and paired landscape inputs were absent from the checked-out repository. |

No external method was redownloaded or rerun. Existing cached predictions were read once by the v1.3 entry point and old v0.x-v1.2 outputs were not overwritten.
"""
    (REPORT / "METHOD_REPRODUCIBILITY.md").write_text(text, encoding="utf-8")


def main():
    ensure_dirs()
    comp = parse_competition()
    audit = write_audit(comp)
    write_contract()
    per = per_protein_decomposition(comp)
    protein_bootstrap_summary(per)
    leave_one_protein_out(per)
    permutation_nulls(comp, per)
    conditionality()
    cached_shuffle_diagnostics()
    figure2(per)
    method_reproducibility()
    update_progress(per, audit)
    print(json.dumps({"status": "complete", "competition_rows": len(comp), "proteins": sorted(comp.protein_id.unique()), "results": OUT.relative_to(ROOT).as_posix()}, indent=2))


if __name__ == "__main__":
    main()
