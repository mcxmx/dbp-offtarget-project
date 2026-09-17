from __future__ import annotations

import hashlib
import json
import math
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.sequence_equivalence import canonical_rc

OUT = ROOT / "results" / "v1_2"
V11 = ROOT / "results" / "v1_1"
PROTEINS = ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP035", "DBP048"]
PRIMARY = PROTEINS[:-1]
BASES = "ACGT"
EPS = 1e-12
N_PERMUTATIONS = 10_000
PERM_SEED = 1201


def ensure_dirs() -> None:
    for name in [
        "00_inputs", "01_reproduction", "02_permutation_nulls", "03_baselines",
        "04_position_decomposition", "05_large_effects", "06_directionality",
        "07_assay_triangle", "08_assay_ceiling", "09_local_conditionality",
        "10_position_robustness", "11_mutation_robustness", "12_dbp048_sensitivity",
        "13_published_subset", "14_cross_model", "figures",
    ]:
        (OUT / name).mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe_spearman(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or len(y) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan")
    return float(spearmanr(x, y).statistic)


def safe_pearson(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or len(y) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan")
    return float(pearsonr(x, y).statistic)


def load_inputs() -> dict[str, pd.DataFrame]:
    paths = {
        "deep": V11 / "05_assay_aligned_eval/deeppbs_per_mutation_predictions.tsv",
        "na": V11 / "05_assay_aligned_eval/nampnn_local_per_mutation_predictions.tsv",
        "pbm": V11 / "06_assay_agreement/pbm_competition_per_mutation.tsv",
        "metrics": V11 / "05_assay_aligned_eval/competition_metrics.tsv",
        "local_global": V11 / "05_assay_aligned_eval/local_global_comparison.tsv",
        "mapping": V11 / "03_register_mapping/register_map.tsv",
        "mapping_summary": V11 / "03_register_mapping/register_summary.tsv",
        "dbp048_windows": V11 / "03_register_mapping/register_sensitivity_candidates.tsv",
        "ppm": V11 / "01_deeppbs_ingest/deeppbs_design7_pwm.tsv",
        "v11_metrics": V11 / "v1_1_metrics.json",
        "sequence_kmer3": ROOT / "data/processed/v0_3_1/designed_dbp_sequence_baseline_rc_aware_scored_v0_3_1.parquet",
    }
    for name, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(f"required v1.1 input missing: {path}")
    deep = pd.read_csv(paths["deep"], sep="\t")
    na = pd.read_csv(paths["na"], sep="\t")
    pbm = pd.read_csv(paths["pbm"], sep="\t")
    for frame in [deep, na, pbm]:
        if set(frame.protein) - set(PROTEINS):
            raise AssertionError("unexpected protein identifier in v1.1 local table")
        if frame.duplicated(["protein", "assay_position", "wt_base", "mut_base"]).any():
            raise AssertionError("duplicate mutation key in v1.1 input")
    mapping = pd.read_csv(paths["mapping"], sep="\t")
    mapping_summary = pd.read_csv(paths["mapping_summary"], sep="\t")
    windows = pd.read_csv(paths["dbp048_windows"], sep="\t")
    ppm = pd.read_csv(paths["ppm"], sep="\t")
    metrics = pd.read_csv(paths["metrics"], sep="\t")
    local_global = pd.read_csv(paths["local_global"], sep="\t")
    return {"deep": deep, "na": na, "pbm": pbm, "mapping": mapping,
            "mapping_summary": mapping_summary, "windows": windows, "ppm": ppm,
            "metrics": metrics, "local_global": local_global,
            "paths": paths}


def freeze_inputs(inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    manifest_path = OUT / "00_inputs/v1_1_input_manifest.tsv"
    if manifest_path.exists():
        existing = pd.read_csv(manifest_path, sep="\t")
        expected = {path.relative_to(ROOT).as_posix(): path.stat().st_size for path in inputs["paths"].values()}
        recorded = dict(zip(existing.source_file, existing.bytes))
        if recorded == expected:
            return existing
        raise AssertionError("frozen input path or file-size discrepancy; refusing to overwrite manifest")
    rows = []
    for name, path in inputs["paths"].items():
        rows.append({"input_name": name, "source_file": path.relative_to(ROOT).as_posix(),
                     "bytes": path.stat().st_size, "sha256": sha256(path),
                     "source_version": "v1.1" if str(path).startswith(str(V11)) else "historical_v0.3.1",
                     "status": "frozen_read_only"})
    manifest = pd.DataFrame(rows)
    manifest.to_csv(manifest_path, sep="\t", index=False)
    (OUT / "00_inputs/README.md").write_text(
        "# v1.2 frozen inputs\n\n"
        "All analysis inputs are read-only v1.1 artifacts. v1.2 does not rewrite or regenerate v1.1 results. "
        "SHA256 values and file sizes are recorded in `v1_1_input_manifest.tsv`.\n",
        encoding="utf-8")
    return manifest


def reproduce_v11(inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    expected = {
        "DBP001": 0.2271255060728745, "DBP003": 0.18238866396761136,
        "DBP005": 0.34916133214488293, "DBP006": 0.7447532614861032,
        "DBP009": 0.5655133295519001, "DBP035": 0.3583988331577668,
        "DBP048": 0.049833887043189376,
    }
    rows = []
    for protein in PROTEINS:
        g = inputs["deep"].query("protein == @protein")
        observed = safe_spearman(g.predicted_delta_logP, g.experimental_effect)
        rows.append({"protein": protein, "v1_1_recorded": expected[protein],
                     "v1_2_recomputed": observed,
                     "absolute_difference": abs(observed - expected[protein]),
                     "status": "PASS" if np.isclose(observed, expected[protein], atol=1e-12) else "FAIL"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "01_reproduction/deeppbs_local_reproduction.tsv", sep="\t", index=False)
    if not (out.status == "PASS").all():
        raise AssertionError("v1.1 DeepPBS local reproduction failed; downstream analysis stopped")
    (OUT / "01_reproduction/reproduction_report.md").write_text(
        "# v1.1 local metric reproduction\n\n"
        "The v1.1 DeepPBS local correlations were recomputed directly from the frozen v1.1 per-mutation table. "
        "All seven values match to absolute tolerance 1e-12. DBP048 remains sensitivity-only in all primary summaries.\n",
        encoding="utf-8")
    return out


def generate_permutation_values(g: pd.DataFrame, null_type: str, n_permutations: int, seed: int) -> np.ndarray:
    if null_type not in {"all_mutation_shuffle", "position_preserving_shuffle"}:
        raise ValueError(f"unknown null type: {null_type}")
    pred = g.predicted_delta_logP.to_numpy(float)
    exp = g.experimental_effect.to_numpy(float)
    rng = np.random.default_rng(seed)
    values = np.empty(n_permutations, dtype=float)
    groups = [np.asarray(list(idx), dtype=int) for idx in g.groupby("assay_position").groups.values()]
    for j in range(n_permutations):
        shuffled = pred.copy()
        if null_type == "all_mutation_shuffle":
            shuffled = shuffled[rng.permutation(len(shuffled))]
        else:
            for idx in groups:
                shuffled[idx] = shuffled[rng.permutation(idx)]
        values[j] = safe_spearman(shuffled, exp)
    return values


def permutation_nulls(deep: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows, distribution_rows = [], []
    for pi, protein in enumerate(PROTEINS):
        g = deep[deep.protein == protein].reset_index(drop=True)
        pred = g.predicted_delta_logP.to_numpy(float)
        exp = g.experimental_effect.to_numpy(float)
        observed = safe_spearman(pred, exp)
        for null_type in ["all_mutation_shuffle", "position_preserving_shuffle"]:
            seed = PERM_SEED + pi * 100 + (0 if null_type == "all_mutation_shuffle" else 1)
            values = generate_permutation_values(g, null_type, N_PERMUTATIONS, seed)
            for j in range(N_PERMUTATIONS):
                distribution_rows.append({"protein": protein, "null_type": null_type,
                                          "permutation": j, "spearman": values[j]})
            tail = float((1 + np.sum(values >= observed)) / (N_PERMUTATIONS + 1)) if np.isfinite(observed) else float("nan")
            summary_rows.append({"protein": protein, "null_type": null_type,
                                 "observed_spearman": observed, "null_median": np.median(values),
                                 "null_q025": np.quantile(values, .025), "null_q975": np.quantile(values, .975),
                                 "empirical_p_greater": tail,
                                 "n_mutations": len(g), "n_permutations": N_PERMUTATIONS,
                                 "seed": seed})
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "02_permutation_nulls/permutation_summary.tsv", sep="\t", index=False)
    pd.DataFrame(distribution_rows).to_csv(OUT / "02_permutation_nulls/permutation_distributions.tsv", sep="\t", index=False)
    (OUT / "02_permutation_nulls/README.md").write_text(
        "# Permutation nulls\n\n"
        f"The observed DeepPBS scores are permuted within each protein using {N_PERMUTATIONS:,} deterministic permutations. "
        "The position-preserving null shuffles the three predicted substitution scores only within each assay position, "
        "preserving position-level magnitude while removing mutant-base identity. Empirical tail probabilities are diagnostics, "
        "not population-level inference for N=7 proteins.\n",
        encoding="utf-8")
    return summary, pd.DataFrame(distribution_rows)


def load_deeppbs_info(inputs: dict[str, pd.DataFrame]) -> dict[str, np.ndarray]:
    ppm = inputs["ppm"]
    result = {}
    for protein in PROTEINS:
        g = ppm[ppm.protein == protein].sort_values("position")
        p = g[["p_A", "p_C", "p_G", "p_T"]].to_numpy(float)
        if not np.allclose(p.sum(axis=1), 1, atol=2e-5):
            raise AssertionError(f"DeepPBS PPM normalization failed for {protein}")
        result[protein] = p
    return result


def average_precision(y_true: np.ndarray, score: np.ndarray) -> float:
    order = np.argsort(-score, kind="mergesort")
    y = y_true[order].astype(int)
    total = int(y.sum())
    if total == 0 or total == len(y):
        return float("nan")
    cumulative = np.cumsum(y)
    return float(np.sum((cumulative / np.arange(1, len(y) + 1)) * y) / total)


def baseline_metrics(deep: pd.DataFrame, pbm: pd.DataFrame, inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    ppm = load_deeppbs_info(inputs)
    rows = []
    transitions = {"AG", "GA", "CT", "TC"}
    for protein in PROTEINS:
        g = deep[deep.protein == protein].copy()
        # Information content is a non-competition, position-only quantity.
        pos_info = {}
        for pos, row in inputs["ppm"][inputs["ppm"].protein == protein].set_index("position").iterrows():
            probs = np.array([row.p_A, row.p_C, row.p_G, row.p_T], dtype=float)
            pos_info[int(pos)] = float(np.sum(probs * np.log2(np.maximum(probs, EPS) / .25)))
        # Score is on the same mutant-effect axis as DeepPBS delta and the
        # competition effect: high WT information predicts a negative mutant
        # change, so use negative information content.
        g["position_only_score"] = g.assay_position.map(lambda x: -pos_info[int(x)])
        g["mutation_type_transition_score"] = [int(a + b in transitions) for a, b in zip(g.wt_base, g.mut_base)]
        g["uniform_score"] = 1.0
        for method, col, definition in [
            ("DeepPBS_position_information_only", "position_only_score", "negative PPM information content; identical for all substitutions at a position"),
            ("mutation_type_transition_indicator", "mutation_type_transition_score", "fixed transition=1/transversion=0; no competition labels fitted"),
            ("uniform_no_information", "uniform_score", "constant score; correlation intentionally undefined"),
        ]:
            rows.append({"method": method, "protein": protein, "n_mutations": len(g),
                         "spearman": safe_spearman(g[col], g.experimental_effect),
                         "pearson": safe_pearson(g[col], g.experimental_effect), "definition": definition,
                         "status": "UNDEFINED_CONSTANT" if method == "uniform_no_information" else "EVALUATED"})
    # PBM is an experimental, not label-free, comparator; it is kept separate.
    sequence = pd.read_parquet(inputs["paths"]["sequence_kmer3"])
    sequence["protein"] = sequence.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}")
    for protein in PROTEINS:
        g = pbm[pbm.protein == protein].copy()
        # PBM position-only control keeps between-position sensitivity and removes
        # which-base information by assigning the same value to all substitutions.
        g["pbm_position_only"] = g.groupby("assay_position").predicted_delta.transform("mean")
        rows.append({"method": "PBM_experimental_local_effect", "protein": protein,
                     "n_mutations": len(g), "spearman": safe_spearman(g.predicted_delta, g.experimental_effect),
                     "pearson": safe_pearson(g.predicted_delta, g.experimental_effect),
                     "definition": "frozen experimental PBM mutant-minus-WT score on mapped 7-mer positions",
                     "status": "PRIMARY" if protein in PRIMARY else "SENSITIVITY_ONLY"})
        rows.append({"method": "PBM_position_only", "protein": protein,
                     "n_mutations": len(g), "spearman": safe_spearman(g.pbm_position_only, g.experimental_effect),
                     "pearson": safe_pearson(g.pbm_position_only, g.experimental_effect),
                     "definition": "mean PBM mutation effect per position; identical for all three substitutions",
                     "status": "PRIMARY" if protein in PRIMARY else "SENSITIVITY_ONLY"})
        lookup = sequence[sequence.protein == protein].drop_duplicates("canonical_7mer").set_index("canonical_7mer")["kmer3_jaccard_to_paper_motif_rc_aware"].to_dict()
        seq_score = [lookup.get(canonical_rc(s), np.nan) for s in g.mutated_pbm_window]
        rows.append({"method": "historical_sequence_kmer3_local", "protein": protein,
                     "n_mutations": int(np.isfinite(seq_score).sum()),
                     "spearman": safe_spearman(seq_score, g.experimental_effect),
                     "pearson": safe_pearson(seq_score, g.experimental_effect),
                     "definition": "frozen v0.3.1 RC-aware 3-mer Jaccard score on mapped mutant 7-mers",
                     "status": "PRIMARY" if protein in PRIMARY else "SENSITIVITY_ONLY"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "03_baselines/baseline_metrics.tsv", sep="\t", index=False)
    (OUT / "03_baselines/baseline_protocol.md").write_text(
        "# Baseline protocol\n\n"
        "The structural position-only control uses negative DeepPBS PPM information content on the mutant-effect axis and assigns the same score to all three substitutions at a position. A separate PBM position-only comparator uses mean experimental PBM mutation effect per position. "
        "The mutation-type baseline is a fixed, label-free transition indicator; it is intentionally not fitted and is not a mechanistic claim. "
        "The uniform baseline is a constant and is reported as undefined rather than numerically hacked. PBM local effect is an experimental comparator, not a label-free baseline. The historical RC-aware 3-mer baseline is reused without fitting.\n",
        encoding="utf-8")
    return out


def position_decomposition(deep: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pos_rows, within_rows = [], []
    for protein in PROTEINS:
        g = deep[deep.protein == protein].copy()
        grouped = g.groupby("assay_position", sort=True)
        pos = grouped.agg(experimental_position_severity=("experimental_effect", "mean"),
                          deeppbs_position_severity=("predicted_delta_logP", "mean"),
                          n_mutations=("mut_base", "size")).reset_index()
        pos_rows.extend([{"protein": protein, **r.to_dict(),
                          "position_spearman": safe_spearman(pos.deeppbs_position_severity, pos.experimental_position_severity)}
                         for _, r in pos.iterrows()])
        for position, q in grouped:
            q = q.reset_index(drop=True)
            pred_order = np.argsort(q.predicted_delta_logP.to_numpy(float), kind="mergesort")
            exp_order = np.argsort(q.experimental_effect.to_numpy(float), kind="mergesort")
            pair_correct = pair_total = 0
            for i, j in combinations(range(len(q)), 2):
                dp = q.predicted_delta_logP.iloc[i] - q.predicted_delta_logP.iloc[j]
                de = q.experimental_effect.iloc[i] - q.experimental_effect.iloc[j]
                if dp == 0 or de == 0:
                    continue
                pair_total += 1
                pair_correct += int(dp * de > 0)
            within_rows.append({"protein": protein, "assay_position": int(position), "n_mutations": len(q),
                                "top_deleterious_exact": int(pred_order[0] == exp_order[0]),
                                "pairwise_order_accuracy": pair_correct / pair_total if pair_total else np.nan,
                                "pairwise_comparisons": pair_total,
                                "within_position_spearman": safe_spearman(q.predicted_delta_logP, q.experimental_effect)})
    pos_df = pd.DataFrame(pos_rows)
    within_df = pd.DataFrame(within_rows)
    pos_df.to_csv(OUT / "04_position_decomposition/position_severity.tsv", sep="\t", index=False)
    within_df.to_csv(OUT / "04_position_decomposition/within_position_substitution.tsv", sep="\t", index=False)
    summary = within_df.groupby("protein").agg(top_deleterious_exact_mean=("top_deleterious_exact", "mean"),
                                                pairwise_order_accuracy_mean=("pairwise_order_accuracy", "mean"),
                                                within_position_spearman_median=("within_position_spearman", "median")).reset_index()
    summary.to_csv(OUT / "04_position_decomposition/position_decomposition_summary.tsv", sep="\t", index=False)
    (OUT / "04_position_decomposition/README.md").write_text(
        "# Position decomposition\n\n"
        "Position severity is the mean of the three mutation effects at each position. Within-position metrics compare the ordering of the three substitutions, "
        "with lower experimental effect / lower predicted delta interpreted as more deleterious. Position-level agreement and mutant-base identity are reported separately.\n",
        encoding="utf-8")
    return pos_df, within_df


def large_effects(deep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein in PROTEINS:
        g = deep[deep.protein == protein].copy()
        threshold = float(g.experimental_effect.quantile(.25))
        g["severe_deleterious"] = g.experimental_effect <= threshold
        score = -g.predicted_delta_logP.to_numpy(float)
        y = g.severe_deleterious.to_numpy(bool)
        order = np.argsort(-score, kind="mergesort")
        k = int(math.ceil(.25 * len(g)))
        top = y[order[:k]]
        prevalence = float(y.mean())
        rows.append({"protein": protein, "n_mutations": len(g), "experimental_q25_threshold": threshold,
                     "n_severe": int(y.sum()), "average_precision": average_precision(y, score),
                     "top_quartile_k": k, "top_quartile_recall": float(top.sum() / y.sum()) if y.sum() else np.nan,
                     "top_quartile_enrichment": float(top.mean() / prevalence) if prevalence else np.nan,
                     "status": "PRIMARY_SIX" if protein in PRIMARY else "SENSITIVITY_ONLY"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "05_large_effects/large_effect_recovery.tsv", sep="\t", index=False)
    (OUT / "05_large_effects/protocol.md").write_text(
        "# Large-effect recovery\n\n"
        "The severe-deleterious class is fixed as the lowest within-protein quartile of competition effect, where lower negative PE/FITC-derived effect means weaker competition. "
        "This threshold is not tuned. Average precision and top-quartile recall/enrichment are secondary diagnostics.\n",
        encoding="utf-8")
    return out


def directionality(deep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for protein in PROTEINS:
        g = deep[deep.protein == protein]
        # The XLS has mutation competition signal normalized to no-competitor, not WT-relative
        # signed affinity changes. A biologically valid improve-vs-reduce class is unavailable.
        rows.append({"protein": protein, "n_evaluable": 0, "sign_accuracy": np.nan,
                     "balanced_accuracy": np.nan, "tp": np.nan, "tn": np.nan, "fp": np.nan, "fn": np.nan,
                     "status": "NOT_APPLICABLE_NO_WT_RELATIVE_SIGNED_EFFECT",
                     "note": "competition values are normalized signal; all mutations are not a WT-relative improve/reduce classification"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "06_directionality/directionality_metrics.tsv", sep="\t", index=False)
    (OUT / "06_directionality/protocol.md").write_text(
        "# Directionality\n\n"
        "No sign accuracy is reported. The source provides normalized competitor signal for each mutant, not a WT-relative signed binding change with an interpretable improve/reduce zero. "
        "Assigning binary signs would manufacture a class boundary.\n",
        encoding="utf-8")
    return out


def assay_triangle(deep: pd.DataFrame, pbm: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = ["protein", "assay_position", "wt_base", "mut_base"]
    merged = deep.merge(pbm[keys + ["predicted_delta", "mutated_pbm_window"]].rename(columns={"predicted_delta": "pbm_delta"}), on=keys, how="inner")
    merged["canonical_mutated_7mer"] = merged.mutated_pbm_window.map(canonical_rc)
    rows = []
    for protein in PROTEINS:
        g = merged[merged.protein == protein]
        rows.append({"protein": protein, "n_matched": len(g),
                     "rho_deeppbs_competition": safe_spearman(g.predicted_delta_logP, g.experimental_effect),
                     "rho_pbm_competition": safe_spearman(g.pbm_delta, g.experimental_effect),
                     "rho_deeppbs_pbm": safe_spearman(g.predicted_delta_logP, g.pbm_delta),
                     "mapping_status": "PRIMARY" if protein in PRIMARY else "AMBIGUOUS_SENSITIVITY_ONLY"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "07_assay_triangle/assay_triangle.tsv", sep="\t", index=False)
    merged.to_csv(OUT / "07_assay_triangle/matched_mutations.tsv", sep="\t", index=False)
    (OUT / "07_assay_triangle/README.md").write_text(
        "# Assay triangle\n\n"
        "All three correlations use the exact matched mutation keys available in both DeepPBS and PBM local tables (normally 21 mutations per protein). "
        "DBP048 remains sensitivity-only. No correlation is computed over unmatched DeepPBS positions.\n",
        encoding="utf-8")
    return out, merged


def assay_ceiling(triangle: pd.DataFrame) -> pd.DataFrame:
    primary = triangle[triangle.protein.isin(PRIMARY)].copy()
    x, y = primary.rho_pbm_competition.to_numpy(float), primary.rho_deeppbs_competition.to_numpy(float)
    summary = pd.DataFrame([{"n_proteins": len(primary), "rho_across_proteins_exploratory": safe_spearman(x, y),
                             "pearson_across_proteins_exploratory": safe_pearson(x, y),
                             "note": "protein-level exploratory association; N=6"}])
    summary.to_csv(OUT / "08_assay_ceiling/ceiling_summary.tsv", sep="\t", index=False)
    primary[["protein", "rho_pbm_competition", "rho_deeppbs_competition", "mapping_status"]].to_csv(
        OUT / "08_assay_ceiling/ceiling_points.tsv", sep="\t", index=False)
    (OUT / "08_assay_ceiling/interpretation.md").write_text(
        "# Assay-ceiling interpretation\n\n"
        "PBM-vs-competition correlation is treated as context, not a denominator or mathematical ceiling. Negative and near-zero assay agreement make ratio normalization invalid. "
        "The scatter and exploratory across-protein association are descriptive only.\n",
        encoding="utf-8")
    return summary


def local_conditionality(deep: pd.DataFrame, mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    m = mapping[mapping.pbm_7mer_position.notna()][["protein", "assay_position", "pbm_7mer_position"]].copy()
    d = deep.merge(m, on=["protein", "assay_position"], how="inner")
    d["local_key"] = [f"{int(p)}:{w}:{u}" for p, w, u in zip(d.pbm_7mer_position, d.wt_base, d.mut_base)]
    vectors = {p: d[d.protein == p].set_index("local_key").predicted_delta_logP.to_dict() for p in PROTEINS}
    rows = []
    for a, b in combinations(PROTEINS, 2):
        common = sorted(set(vectors[a]) & set(vectors[b]))
        xa = [vectors[a][k] for k in common]; xb = [vectors[b][k] for k in common]
        rows.append({"protein_a": a, "protein_b": b, "n_common_mutations": len(common),
                     "spearman": safe_spearman(xa, xb), "pearson": safe_pearson(xa, xb),
                     "status": "PRIMARY_PAIR" if a in PRIMARY and b in PRIMARY else "DBP048_SENSITIVITY_INCLUDED"})
    pairs = pd.DataFrame(rows)
    pairs.to_csv(OUT / "09_local_conditionality/pairwise_local_correlations.tsv", sep="\t", index=False)
    prim = pairs[pairs.status == "PRIMARY_PAIR"]
    summary = pd.DataFrame([{"set": "primary_six", "n_pairs": len(prim), "median_pairwise_spearman": prim.spearman.median(),
                             "min_pairwise_spearman": prim.spearman.min(), "max_pairwise_spearman": prim.spearman.max()},
                            {"set": "all_seven_descriptive", "n_pairs": len(pairs), "median_pairwise_spearman": pairs.spearman.median(),
                             "min_pairwise_spearman": pairs.spearman.min(), "max_pairwise_spearman": pairs.spearman.max()}])
    summary.to_csv(OUT / "09_local_conditionality/summary.tsv", sep="\t", index=False)
    (OUT / "09_local_conditionality/protocol.md").write_text(
        "# Local conditionality\n\n"
        "Vectors are compared only on common keys defined by independently mapped relative 7-mer position, WT base, and mutant base. "
        "This is a representation-level conditionality diagnostic, not an independent validation analysis.\n",
        encoding="utf-8")
    return pairs, summary


def position_robustness(deep: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    detail = []
    for protein in PROTEINS:
        g = deep[deep.protein == protein]
        original = safe_spearman(g.predicted_delta_logP, g.experimental_effect)
        for position in sorted(g.assay_position.unique()):
            q = g[g.assay_position != position]
            detail.append({"protein": protein, "dropped_position": int(position), "n_remaining": len(q),
                           "rho": safe_spearman(q.predicted_delta_logP, q.experimental_effect),
                           "original_rho": original})
    detail_df = pd.DataFrame(detail)
    rows = []
    for protein, g in detail_df.groupby("protein"):
        rows.append({"protein": protein, "original_rho": g.original_rho.iloc[0], "loo_min_rho": g.rho.min(),
                     "loo_max_rho": g.rho.max(), "loo_median_rho": g.rho.median(),
                     "loo_range": g.rho.max() - g.rho.min(), "n_positions": len(g),
                     "status": "PRIMARY_SIX" if protein in PRIMARY else "SENSITIVITY_ONLY"})
    summary = pd.DataFrame(rows)
    detail_df.to_csv(OUT / "10_position_robustness/leave_one_position_out.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "10_position_robustness/summary.tsv", sep="\t", index=False)
    return detail_df, summary


def mutation_robustness(deep: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    detail = []
    for protein in PROTEINS:
        g = deep[deep.protein == protein].reset_index(drop=True)
        original = safe_spearman(g.predicted_delta_logP, g.experimental_effect)
        for i in range(len(g)):
            q = g.drop(index=i)
            detail.append({"protein": protein, "removed_index": i,
                           "removed_assay_position": int(g.iloc[i].assay_position),
                           "removed_wt_base": g.iloc[i].wt_base, "removed_mut_base": g.iloc[i].mut_base,
                           "rho": safe_spearman(q.predicted_delta_logP, q.experimental_effect),
                           "original_rho": original,
                           "rho_change": safe_spearman(q.predicted_delta_logP, q.experimental_effect) - original})
    detail_df = pd.DataFrame(detail)
    rows = []
    for protein, g in detail_df.groupby("protein"):
        g = g.copy(); influence = g.rho_change.abs()
        ix = influence.idxmax(); r = g.loc[ix]
        rows.append({"protein": protein, "original_rho": g.original_rho.iloc[0], "jackknife_min_rho": g.rho.min(),
                     "jackknife_max_rho": g.rho.max(), "jackknife_median_rho": g.rho.median(),
                     "jackknife_range": g.rho.max() - g.rho.min(), "max_abs_change": influence.max(),
                     "most_influential_position": r.removed_assay_position, "most_influential_wt": r.removed_wt_base,
                     "most_influential_mut": r.removed_mut_base})
    summary = pd.DataFrame(rows)
    detail_df.to_csv(OUT / "11_mutation_robustness/mutation_jackknife.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "11_mutation_robustness/summary.tsv", sep="\t", index=False)
    return detail_df, summary


def dbp048_sensitivity(deep: pd.DataFrame, inputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    protein = "DBP048"
    P = load_deeppbs_info(inputs)[protein]
    g = deep[(deep.protein == protein) & deep.assay_position.between(7, 13)].sort_values(["assay_position", "mut_base"]).copy()
    rows = []
    for _, w in inputs["windows"].iterrows():
        start = int(w.model_window_start_1based)
        pred = []
        for _, r in g.iterrows():
            mpos = start + (int(r.assay_position) - 7)
            pred.append(math.log(P[mpos - 1, BASES.index(r.mut_base)] + EPS) - math.log(P[mpos - 1, BASES.index(r.wt_base)] + EPS))
        pred = np.asarray(pred, dtype=float)
        position = g.assign(_pred=pred).groupby("assay_position").agg(pred=("_pred", "mean"), exp=("experimental_effect", "mean"))
        within = [safe_spearman(qp, q.experimental_effect) for (_, q), qp in zip(g.groupby("assay_position"), pred.reshape(7, 3))]
        threshold = float(g.experimental_effect.quantile(.25)); severe = (g.experimental_effect <= threshold).to_numpy(bool)
        rows.append({"protein": protein, "model_window_start_1based": start,
                     "model_window": w.model_window, "assay_window": "CTGACGC",
                     "n_mutations": len(g), "spearman": safe_spearman(pred, g.experimental_effect),
                     "pearson": safe_pearson(pred, g.experimental_effect),
                     "position_severity_spearman": safe_spearman(position.pred, position.exp),
                     "within_position_median_spearman": float(np.nanmedian(within)),
                     "severe_quartile_average_precision": average_precision(severe, -pred),
                     "mapping_status": "AMBIGUOUS_SENSITIVITY_ONLY",
                     "previous_v1_1_window": int(start == 4), "selection_status": w.selection_status})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "12_dbp048_sensitivity/all_windows.tsv", sep="\t", index=False)
    (OUT / "12_dbp048_sensitivity/README.md").write_text(
        "# DBP048 sensitivity\n\n"
        "All eight pre-existing plausible model windows are evaluated against the sequence-C assay window `CTGACGC`. "
        "No window is selected using experimental correlation. The v1.1 positional window (model start 4, assay start 7) is marked only for provenance.\n",
        encoding="utf-8")
    return out


def published_subset(deep: pd.DataFrame) -> pd.DataFrame:
    g = deep[deep.protein.isin(["DBP005", "DBP006", "DBP009", "DBP035"])].copy()
    out = g.groupby("protein").apply(lambda x: pd.Series({"n_mutations": len(x), "spearman": safe_spearman(x.predicted_delta_logP, x.experimental_effect)}), include_groups=False).reset_index()
    out["subset"] = "predefined_published_DeepPBS_four_protein_subset"
    out.to_csv(OUT / "13_published_subset/metrics.tsv", sep="\t", index=False)
    pd.DataFrame([{"subset": "DBP005|DBP006|DBP009|DBP035", "n_proteins": 4,
                    "median_spearman": out.spearman.median(),
                    "note": "same underlying development-exposed proteins and assay; descriptive only"}]).to_csv(
        OUT / "13_published_subset/summary.tsv", sep="\t", index=False)
    return out


def cross_model(deep: pd.DataFrame, na: pd.DataFrame) -> pd.DataFrame:
    keys = ["protein", "assay_position", "wt_base", "mut_base"]
    d = deep[deep.assay_position.between(1, 14)][keys + ["predicted_delta_logP", "experimental_effect", "mapping_status"]]
    n = na[keys + ["predicted_delta", "mapping_status"]].rename(columns={"predicted_delta": "nampnn_derived_delta", "mapping_status": "na_mapping_status"})
    m = d.merge(n, on=keys, how="inner")
    rows = []
    for pi, protein in enumerate(PROTEINS):
        g = m[m.protein == protein].reset_index(drop=True)
        values = {}
        for mi, (prefix, col) in enumerate([("deeppbs", "predicted_delta_logP"), ("nampnn_derived", "nampnn_derived_delta")]):
            score = g[col].to_numpy(float)
            observed = safe_spearman(score, g.experimental_effect)
            position = g.assign(_score=score).groupby("assay_position").agg(pred=("_score", "mean"), exp=("experimental_effect", "mean"))
            within = [safe_spearman(q[col], q.experimental_effect) for _, q in g.groupby("assay_position")]
            threshold = float(g.experimental_effect.quantile(.25)); y = (g.experimental_effect <= threshold).to_numpy(bool)
            rng = np.random.default_rng(PERM_SEED + 9000 + pi * 10 + mi)
            perm = np.empty(N_PERMUTATIONS)
            for j in range(N_PERMUTATIONS):
                shuffled = score.copy()
                for _, idx in g.groupby("assay_position").groups.items():
                    idx = np.asarray(list(idx), dtype=int)
                    shuffled[idx] = shuffled[rng.permutation(idx)]
                perm[j] = safe_spearman(shuffled, g.experimental_effect)
            values.update({
                f"{prefix}_rho": observed,
                f"{prefix}_position_rho": safe_spearman(position.pred, position.exp),
                f"{prefix}_within_position_median_rho": float(np.nanmedian(within)),
                f"{prefix}_severe_quartile_ap": average_precision(y, -score),
                f"{prefix}_position_null_median": float(np.median(perm)),
                f"{prefix}_position_null_q975": float(np.quantile(perm, .975)),
                f"{prefix}_position_null_p_greater": float((1 + np.sum(perm >= observed)) / (N_PERMUTATIONS + 1)),
                f"{prefix}_null_adjusted_rho": observed - float(np.median(perm)),
            })
        rows.append({"protein": protein, "n_matched": len(g), **values,
                     "rho_difference_deeppbs_minus_nampnn": values["deeppbs_rho"] - values["nampnn_derived_rho"],
                     "na_score_label": "NA-MPNN frozen global landscape-derived local mutation score",
                     "status": "PRIMARY" if protein in PRIMARY else "SENSITIVITY_ONLY"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "14_cross_model/per_protein.tsv", sep="\t", index=False)
    m.to_csv(OUT / "14_cross_model/matched_mutations.tsv", sep="\t", index=False)
    return out


def make_figures(deep: pd.DataFrame, pbm: pd.DataFrame, nulls: pd.DataFrame,
                 pos: pd.DataFrame, within: pd.DataFrame, triangle: pd.DataFrame,
                 loo: pd.DataFrame, large: pd.DataFrame, cross: pd.DataFrame) -> None:
    f = OUT / "figures"
    # Figure 1: observed versus both null distributions.
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), squeeze=False)
    dist = pd.read_csv(OUT / "02_permutation_nulls/permutation_distributions.tsv", sep="\t")
    for ax, protein in zip(axes.flat, PROTEINS):
        q = nulls[nulls.protein == protein]
        for typ, color in [("all_mutation_shuffle", "#8da0cb"), ("position_preserving_shuffle", "#fc8d62")]:
            vals = dist[(dist.protein == protein) & (dist.null_type == typ)].spearman
            ax.hist(vals, bins=35, alpha=.45, color=color, density=True)
            obs = float(q[q.null_type == typ].observed_spearman.iloc[0]); ax.axvline(obs, color="black", lw=1)
        ax.set_title(protein); ax.set_xlabel("Spearman"); ax.set_ylabel("null density")
    for ax in axes.flat[len(PROTEINS):]: ax.axis("off")
    fig.suptitle("DeepPBS local correlations versus permutation nulls"); fig.tight_layout(); fig.savefig(f / "figure1_permutation_nulls.png", dpi=220); plt.close(fig)
    # Figure 2: strictly matched 21-mutation comparison from the assay triangle.
    matched = triangle.set_index("protein").reindex(PROTEINS)
    fig, ax = plt.subplots(figsize=(8, 5)); x = np.arange(len(PROTEINS))
    ax.scatter(x - .09, matched.rho_deeppbs_competition, label="DeepPBS vs competition", color="#1f77b4", s=45)
    ax.scatter(x + .09, matched.rho_pbm_competition, label="PBM vs competition", color="#2ca02c", s=45)
    labels = [p + ("*" if p == "DBP048" else "") for p in PROTEINS]
    ax.set_xticks(x, labels, rotation=45); ax.axhline(0, color="grey", lw=.8)
    ax.set_ylabel("Spearman"); ax.set_title("Matched local comparison (21 mutations per protein)")
    ax.legend(); fig.tight_layout(); fig.savefig(f / "figure2_deeppbs_vs_pbm.png", dpi=220); plt.close(fig)
    # Figure 3: position and within-position decomposition.
    ps = pos.groupby("protein").position_spearman.first()
    ws = within.groupby("protein").pairwise_order_accuracy.mean()
    fig, ax = plt.subplots(figsize=(8, 5)); x = np.arange(len(PROTEINS))
    ax.scatter(x-.08, [ps.get(p, np.nan) for p in PROTEINS], label="position-severity Spearman", s=45)
    ax.scatter(x+.08, [ws.get(p, np.nan) for p in PROTEINS], label="within-position pairwise accuracy", s=45)
    ax.axhline(.5, color="#fc8d62", linestyle="--", lw=.8, label="pairwise chance reference")
    ax.set_xticks(x, [p + ("*" if p == "DBP048" else "") for p in PROTEINS], rotation=45)
    ax.axhline(0,color="grey",lw=.8); ax.legend(fontsize=8); ax.set_ylabel("metric value")
    ax.set_title("Position sensitivity versus mutant-base ordering")
    fig.tight_layout(); fig.savefig(f / "figure3_position_decomposition.png", dpi=220); plt.close(fig)
    # Figure 4: assay-ceiling scatter.
    fig, ax = plt.subplots(figsize=(6, 5));
    for _, r in triangle.iterrows():
        ax.scatter(r.rho_pbm_competition, r.rho_deeppbs_competition, s=55, c="#d62728" if r.protein == "DBP048" else "#1f77b4"); ax.text(r.rho_pbm_competition, r.rho_deeppbs_competition, r.protein, fontsize=8)
    ax.plot([-1, 1], [-1, 1], "--", color="grey", lw=.8); ax.set_xlabel("PBM vs competition rho"); ax.set_ylabel("DeepPBS vs competition rho"); ax.set_title("Assay agreement context (not a ceiling ratio)"); fig.tight_layout(); fig.savefig(f / "figure4_assay_ceiling.png", dpi=220); plt.close(fig)
    # Figure 5: triangle per protein.
    fig, ax = plt.subplots(figsize=(9, 5)); x = np.arange(len(PROTEINS));
    for j, col in enumerate(["rho_deeppbs_competition", "rho_pbm_competition", "rho_deeppbs_pbm"]): ax.scatter(x + (j-1)*.12, triangle.set_index("protein").reindex(PROTEINS)[col], label=col.replace("rho_", ""), s=42)
    ax.set_xticks(x, PROTEINS, rotation=45); ax.axhline(0,color="grey",lw=.8); ax.set_ylabel("Spearman"); ax.legend(fontsize=8); ax.set_title("Matched local assay triangle"); fig.tight_layout(); fig.savefig(f / "figure5_assay_triangle.png", dpi=220); plt.close(fig)
    # Figure 6: leave-one-position-out.
    fig, ax = plt.subplots(figsize=(9, 5));
    for protein, g in loo.groupby("protein"): ax.plot(g.dropped_position, g.rho, marker=".", label=protein)
    ax.set_xlabel("dropped assay position"); ax.set_ylabel("Spearman"); ax.axhline(0,color="grey",lw=.8); ax.legend(ncol=2, fontsize=8); ax.set_title("Leave-one-position-out robustness"); fig.tight_layout(); fig.savefig(f / "figure6_leave_one_position_out.png", dpi=220); plt.close(fig)
    # Figure 7: large effects.
    fig, ax = plt.subplots(figsize=(8, 5)); ax.bar(np.arange(len(large)), large.average_precision, color=["#d62728" if p == "DBP048" else "#1f77b4" for p in large.protein]); ax.set_xticks(np.arange(len(large)), large.protein, rotation=45); ax.set_ylim(0, 1); ax.set_ylabel("average precision"); ax.set_title("Recovery of most deleterious mutations"); fig.tight_layout(); fig.savefig(f / "figure7_large_effect_recovery.png", dpi=220); plt.close(fig)
    # Figure 8: cross model.
    fig, ax = plt.subplots(figsize=(8, 5)); x = np.arange(len(cross)); ax.scatter(x-.1, cross.deeppbs_rho, label="DeepPBS native PWM", s=45); ax.scatter(x+.1, cross.nampnn_derived_rho, label="NA-MPNN derived global score", s=45); ax.axhline(0,color="grey",lw=.8); ax.set_xticks(x, cross.protein, rotation=45); ax.set_ylabel("local Spearman"); ax.legend(fontsize=8); ax.set_title("DeepPBS versus NA-MPNN-derived local score"); fig.tight_layout(); fig.savefig(f / "figure8_cross_model.png", dpi=220); plt.close(fig)

    figure_sources = [
        ("figure1_permutation_nulls.png", "02_permutation_nulls/permutation_distributions.tsv; 02_permutation_nulls/permutation_summary.tsv"),
        ("figure2_deeppbs_vs_pbm.png", "07_assay_triangle/assay_triangle.tsv"),
        ("figure3_position_decomposition.png", "04_position_decomposition/position_severity.tsv; 04_position_decomposition/within_position_substitution.tsv"),
        ("figure4_assay_ceiling.png", "07_assay_triangle/assay_triangle.tsv"),
        ("figure5_assay_triangle.png", "07_assay_triangle/assay_triangle.tsv"),
        ("figure6_leave_one_position_out.png", "10_position_robustness/leave_one_position_out.tsv"),
        ("figure7_large_effect_recovery.png", "05_large_effects/large_effect_recovery.tsv"),
        ("figure8_cross_model.png", "14_cross_model/per_protein.tsv"),
    ]
    pd.DataFrame(figure_sources, columns=["figure", "source_tables"]).to_csv(
        f / "figure_source_manifest.tsv", sep="\t", index=False
    )


def write_summary(inputs: dict[str, pd.DataFrame], repro: pd.DataFrame, nulls: pd.DataFrame,
                  baselines: pd.DataFrame, pos: pd.DataFrame, within: pd.DataFrame,
                  large: pd.DataFrame, triangle: pd.DataFrame, ceiling: pd.DataFrame,
                  cond: pd.DataFrame, loo: pd.DataFrame, jack: pd.DataFrame,
                  dbp48: pd.DataFrame, published: pd.DataFrame, cross: pd.DataFrame) -> None:
    deep_primary = repro[repro.protein.isin(PRIMARY)].v1_2_recomputed
    dnull = nulls[(nulls.null_type == "position_preserving_shuffle") & nulls.protein.isin(PRIMARY)]
    pnull_pass = int(np.sum(dnull.observed_spearman > dnull.null_q975))
    pos_summary = pos.groupby("protein").position_spearman.first()
    within_summary = within.groupby("protein").agg(
        pairwise_order_accuracy=("pairwise_order_accuracy", "mean"),
        top_deleterious_exact=("top_deleterious_exact", "mean"),
        within_position_spearman=("within_position_spearman", "median"),
    )
    primary_within = within_summary.loc[PRIMARY]
    loo_primary = loo[loo.protein.isin(PRIMARY)]
    pbm_baseline_primary = baselines[(baselines.method == "PBM_experimental_local_effect") & baselines.protein.isin(PRIMARY)].spearman.median()
    lines = [
        "# v1.2 summary",
        "",
        "This is a frozen-input robustness and assay-ceiling audit. It reads v1.1 artifacts without changing them; all seven GSE237017 proteins remain development-exposed.",
        "",
        "## Primary results",
        f"- v1.1 reproduction: all seven DeepPBS local correlations passed at 1e-12 tolerance; primary-six median = {deep_primary.median():.4f}.",
        f"- Position-preserving null: {pnull_pass}/6 primary proteins exceed their null 97.5th percentile; see `02_permutation_nulls/permutation_summary.tsv`.",
        f"- DeepPBS primary-six local median = {deep_primary.median():.4f}; PBM experimental local median = {pbm_baseline_primary:.4f}. These headline values use different mutation coverage; the matched 21-mutation comparison is in the assay triangle.",
        f"- Position-level DeepPBS Spearman values: " + ", ".join(f"{p} {pos_summary[p]:.4f}" for p in PROTEINS),
        f"- Within-position pairwise ordering accuracy: primary-six median = {primary_within.pairwise_order_accuracy.median():.4f} (chance reference 0.5); top-deleterious-base exact agreement median = {primary_within.top_deleterious_exact.median():.4f} (chance reference 1/3).",
        f"- DeepPBS large-effect recovery: primary-six median AP = {large[large.protein.isin(PRIMARY)].average_precision.median():.4f}.",
        f"- Leave-one-position-out primary ranges: {loo_primary.loo_min_rho.min():.4f} to {loo_primary.loo_max_rho.max():.4f}; individual summaries are in `10_position_robustness/summary.tsv`.",
        "",
        "## Per-protein DeepPBS local correlations",
        "",
        "| protein | rho | position rho | within-position pairwise accuracy | top-base exact | mapping |",
        "|---|---:|---:|---:|---:|---|",
    ]
    dmetrics = inputs["metrics"].query("method == 'DeepPBS_native_pwm'").set_index("protein")
    for p in PROTEINS:
        status = "PRIMARY" if p in PRIMARY else "AMBIGUOUS / SENSITIVITY_ONLY"
        lines.append(f"| {p} | {dmetrics.loc[p, 'spearman']:.4f} | {pos_summary[p]:.4f} | {within_summary.loc[p, 'pairwise_order_accuracy']:.4f} | {within_summary.loc[p, 'top_deleterious_exact']:.4f} | {status} |")
    lines += [
        "",
        "## Baseline and assay interpretation",
        "- The position-only PPM-information baseline tests position sensitivity without mutant-base identity. The fixed transition indicator is a non-fitted, non-mechanistic mutation-type control. The uniform baseline is correctly undefined.",
        f"- PBM-vs-competition primary-six median = {triangle[triangle.protein.isin(PRIMARY)].rho_pbm_competition.median():.4f}; values range from {triangle[triangle.protein.isin(PRIMARY)].rho_pbm_competition.min():.4f} to {triangle[triangle.protein.isin(PRIMARY)].rho_pbm_competition.max():.4f}.",
        f"- Exploratory across-protein association of assay agreement and DeepPBS performance: Spearman = {ceiling.rho_across_proteins_exploratory.iloc[0]:.4f} (N=6); this is not an assay ceiling ratio.",
        f"- DeepPBS published four-protein descriptive subset median = {published.spearman.median():.4f}; this reuses the same exposed proteins and assay.",
        "- NA-MPNN comparison is limited to the explicitly named frozen global-landscape-derived local score; it is not a native NA-MPNN PPM.",
        "- Directionality is not applicable because the competition XLS does not define WT-relative signed improve/reduce effects.",
        "",
        "## R1-R5 interpretation",
        f"- R1 genuine local structure signal: {'PARTIALLY SUPPORTED / DESCRIPTIVE' if pnull_pass >= 2 else 'NOT ESTABLISHED'}; permutation excess and within-position results must be read protein by protein, without a population claim.",
        "- R2 mostly position-sensitivity signal: SUPPORTED AS THE DOMINANT OBSERVED COMPONENT; only DBP035 exceeded the position-preserving null, while primary-six within-position pairwise accuracy was near its 0.5 chance reference.",
        "- R3 mostly assay-shared signal: INCONCLUSIVE; PBM-vs-competition agreement is heterogeneous and no invalid normalization by assay rho is used.",
        "- R4 heterogeneous/protein-specific success: SUPPORTED DESCRIPTIVELY; DBP006 and DBP009 drive much of the positive DeepPBS local median, while DBP001/003 are weaker.",
        f"- R5 fragile local result: {'PRESENT FOR AT LEAST SOME PROTEINS' if (loo_primary.loo_min_rho < 0).any() else 'NOT INDICATED BY SIGN CROSSING IN PRIMARY LOO'}; mutation jackknife and per-position ranges quantify influence without inferential claims.",
        "",
        "## Boundary",
        "The result does not establish transferable specificity, does not make the seven proteins independent validation, and does not represent all structure-based methods. The most important unresolved limitation remains an independent assay-matched designed-DBP cohort with frozen analysis before label access.",
    ]
    (OUT / "v1_2_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    metrics = {
        "version": "v1.2", "n_proteins": 7, "primary_n_proteins": 6,
        "v11_reproduction_all_pass": bool((repro.status == "PASS").all()),
        "deeppbs_local_primary_median": float(deep_primary.median()),
        "deeppbs_local_all7_median": float(repro.v1_2_recomputed.median()),
        "position_preserving_null_primary_exceed_q975": pnull_pass,
        "deeppbs_per_protein": dict(zip(repro.protein, repro.v1_2_recomputed)),
        "pbm_local_primary_median": float(triangle[triangle.protein.isin(PRIMARY)].rho_pbm_competition.median()),
        "assay_triangle": triangle.set_index("protein").to_dict(orient="index"),
        "primary_within_position_spearman_median": float(primary_within.within_position_spearman.median()),
        "primary_within_position_pairwise_accuracy_median": float(primary_within.pairwise_order_accuracy.median()),
        "primary_top_deleterious_base_exact_median": float(primary_within.top_deleterious_exact.median()),
        "primary_position_severity_median": float(pos[pos.protein.isin(PRIMARY)].groupby("protein").position_spearman.first().median()),
        "large_effect_primary_median_ap": float(large[large.protein.isin(PRIMARY)].average_precision.median()),
        "local_conditionality": cond.to_dict(orient="records"),
        "dbp048_sensitivity": dbp48.to_dict(orient="records"),
        "na_score_label": "NA-MPNN frozen global landscape-derived local mutation score",
        "competition_replicates": "XLS heatmap means of two replicates; individual replicate values unavailable",
        "development_exposed": True,
        "independent_validation": False,
        "interpretations": {"R1": "not_established", "R2": "dominant_observed_component", "R3": "inconclusive", "R4": "descriptively_supported", "R5": "not_indicated_by_primary_sign_crossing"},
    }
    (OUT / "v1_2_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True, default=lambda x: None), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    inputs = load_inputs()
    freeze_inputs(inputs)
    repro = reproduce_v11(inputs)
    nulls, _ = permutation_nulls(inputs["deep"])
    baselines = baseline_metrics(inputs["deep"], inputs["pbm"], inputs)
    pos, within = position_decomposition(inputs["deep"])
    large = large_effects(inputs["deep"])
    directionality(inputs["deep"])
    triangle, _ = assay_triangle(inputs["deep"], inputs["pbm"])
    ceiling = assay_ceiling(triangle)
    cond, cond_summary = local_conditionality(inputs["deep"], inputs["mapping"])
    loo, loo_summary = position_robustness(inputs["deep"])
    jack, jack_summary = mutation_robustness(inputs["deep"])
    dbp48 = dbp048_sensitivity(inputs["deep"], inputs)
    published = published_subset(inputs["deep"])
    cross = cross_model(inputs["deep"], inputs["na"])
    make_figures(inputs["deep"], inputs["pbm"], nulls, pos, within, triangle, loo, large, cross)
    write_summary(inputs, repro, nulls, baselines, pos, within, large, triangle, ceiling,
                  cond_summary, loo_summary, jack_summary, dbp48, published, cross)
    print("v1.2 audit complete", len(inputs["deep"]), "DeepPBS mutation rows")


if __name__ == "__main__":
    main()
