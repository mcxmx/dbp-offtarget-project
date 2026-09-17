"""v0.7 generalization-gap benchmark and failure-analysis artifacts.

This module only reads frozen v0.3-v0.6 artifacts. It does not train a model,
read external quantitative labels, or modify historical result directories.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from difflib import SequenceMatcher

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0_7_benchmark"
DOC = ROOT / "docs" / "v0_7_benchmark"
EXT = ROOT / "metadata" / "v0_7_external"
META = ROOT / "metadata" / "v0_7_benchmark"
SEEDS = "17|29|43"
DBPS = ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]


def safe_spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(spearmanr(a, b).statistic)


def fmt(v, n=4):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "NA"
    if isinstance(v, (np.integer, int)):
        return str(int(v))
    if isinstance(v, (np.floating, float)):
        return f"{float(v):.{n}f}"
    return str(v)


def md_table(df):
    lines = ["| " + " | ".join(df.columns) + " |", "| " + " | ".join("---" for _ in df.columns) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(fmt(x) for x in row) + " |")
    return "\n".join(lines)


def write(df, name):
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / name, index=False)


def per_protein_string(df, id_col="dbp_id", value_col="spearman"):
    if df.empty:
        return ""
    return "|".join(f"{r[id_col]}:{float(r[value_col]):.6f}" for _, r in df.iterrows())


def build_master():
    rows = []
    # Matched v0.5 primary and strict results.
    p05 = pd.read_csv(ROOT / "results/v0_5/primary_per_protein_results.csv")
    m05 = pd.read_csv(ROOT / "results/v0_5/primary_macro_summary.csv")
    for protocol, macro_file, per_file, exposure in [
        ("protein_cluster_loco", "results/v0_5/primary_macro_summary.csv", "results/v0_5/primary_per_protein_results.csv", "development_exposed_all7"),
        ("combined_component_loco", "results/v0_5/strict_component_macro_summary.csv", "results/v0_5/strict_component_per_protein_results.csv", "development_exposed_all7"),
    ]:
        macro = pd.read_csv(ROOT / macro_file)
        per = pd.read_csv(ROOT / per_file)
        for _, r in macro.iterrows():
            vals = per[r["model"]].to_numpy(dtype=float)
            rows.append({
                "stage": "v0.5", "method": r["model"], "model_family": "conditional_pairwise_head",
                "dataset": "GSE237017 designed DBP", "training_proteins": "protein clusters excluding held-out fold",
                "held_out_proteins": "7 proteins (fold-specific)", "assay": "uPBM processed PBM E-score",
                "dna_kmer_length": 7, "split_protocol": protocol, "exposure_status": exposure,
                "model_selection_exposure": "frozen_primary_no_test_selection", "median_spearman": r.get("all7_macro_median", r.get("all_evaluated_macro_median")),
                "per_protein_spearman": per_protein_string(per[["dbp_id", r["model"]]].rename(columns={r["model"]: "spearman"}), value_col="spearman"),
                "random_seeds": r["seeds"], "protein_conditionality_pass": False, "target_conditionality_pass": False,
                "evaluation_role": "developmental benchmark", "comparability_status": "matched within stage",
                "notes": "v0.5 frozen; strict rows use combined-component sensitivity protocol"})
    # v0.6 corrected replay.
    p06 = pd.read_csv(ROOT / "results/v0_6_representation/v06_per_protein_spearman.csv")
    for model in ["M0", "M1", "M1c", "M2", "M3"]:
        vals = p06[model]
        rows.append({
            "stage": "v0.6", "method": model, "model_family": "corrected_rc_ordered_target_pairwise_head",
            "dataset": "GSE237017 designed DBP", "training_proteins": "protein clusters excluding held-out fold",
            "held_out_proteins": "7 proteins (fold-specific)", "assay": "uPBM processed PBM E-score", "dna_kmer_length": 7,
            "split_protocol": "protein_cluster_loco", "exposure_status": "development_exposed_all7",
            "model_selection_exposure": "fixed_replay_no_search", "median_spearman": float(vals.median()),
            "per_protein_spearman": per_protein_string(p06[["dbp_id", model]].rename(columns={model: "spearman"}), value_col="spearman"),
            "random_seeds": SEEDS, "protein_conditionality_pass": False, "target_conditionality_pass": False,
            "evaluation_role": "developmental representation replay", "comparability_status": "matched to v0.5 objective/split",
            "notes": "RC-class injective candidate and ordered target windows"})
    # Historical sequence-only baseline.
    seq = pd.read_csv(ROOT / "results/v0_3_1/tables/designed_dbp_sequence_baseline_rc_aware.csv")
    seq = seq[seq.metric == "kmer3_jaccard_to_paper_motif_rc_aware"]
    rows.append({
        "stage": "v0.3.1", "method": "sequence_kmer3", "model_family": "sequence-only k-mer similarity",
        "dataset": "GSE237017 designed DBP", "training_proteins": "none", "held_out_proteins": "7", "assay": "uPBM processed PBM E-score",
        "dna_kmer_length": 7, "split_protocol": "per-protein ranking/no training", "exposure_status": "development_exposed_all7",
        "model_selection_exposure": "historical_fixed", "median_spearman": float(seq.spearman.median()),
        "per_protein_spearman": per_protein_string(seq.rename(columns={"protein_id": "dbp_id"}), id_col="dbp_id", value_col="spearman"),
        "random_seeds": "none", "protein_conditionality_pass": False, "target_conditionality_pass": False,
        "evaluation_role": "sequence-only comparator", "comparability_status": "historical comparator; no learned protein conditioning",
        "notes": "RC-aware k-mer similarity; not an independent validation"})
    # Natural SimplePC context is deliberately marked unmatched.
    spc = pd.read_csv(ROOT / "results/v0_4_1/tables/simple_pc_performance.csv")
    for p, f, dataset, med, prot, role in [
        ("SimpleProteinConditionalBaseline", "composition_ridge", "natural PBM UniPROBE", 0.3013536821293048, "9 natural proteins", "natural held-out context"),
        ("SimpleProteinConditionalBaseline", "composition_ridge", "GSE237017 designed DBP", 0.3616108129479222, "7 designed proteins", "natural-to-designed context"),
    ]:
        subset = spc[spc.dataset.eq("natural_test" if "natural" in dataset else "designed_external")]
        rows.append({"stage": "v0.4.1", "method": p, "model_family": f, "dataset": dataset,
                     "training_proteins": "natural PBM train split", "held_out_proteins": prot,
                     "assay": "UniPROBE contiguous 8-mer E-score" if "natural" in dataset else "uPBM processed PBM E-score",
                     "dna_kmer_length": 8 if "natural" in dataset else 7, "split_protocol": "v0.4.1 natural 40% cluster split; unmatched external designed evaluation",
                     "exposure_status": "development_exposed_historical", "model_selection_exposure": "historical_v0.4.1",
                     "median_spearman": med, "per_protein_spearman": per_protein_string(subset, id_col="protein_id", value_col="spearman"),
                     "random_seeds": "not reported", "protein_conditionality_pass": "not assessed", "target_conditionality_pass": "not assessed",
                     "evaluation_role": role, "comparability_status": "NOT_MATCHED_PROTOCOL", "notes": "Cannot be used as a causal transfer result or independent validation"})
    # Local/residue smoke.
    local = pd.read_csv(ROOT / "results/v0_5_local/smoke_per_protein_results.csv")
    shuf_local = pd.read_csv(ROOT / "results/v0_5_local/shuffle_diagnostics.csv")
    for model, g in local.groupby("model"):
        test = g[g.partition == "test"]
        pc = shuf_local[shuf_local.model == model]
        rows.append({"stage": "v0.5-local", "method": model, "model_family": "residue/local interaction smoke",
                     "dataset": "GSE237017 designed DBP", "training_proteins": "fold 2 training proteins", "held_out_proteins": "DBP35|DBP5",
                     "assay": "uPBM processed PBM E-score", "dna_kmer_length": 7, "split_protocol": "one-fold protein_cluster_loco_fold_2",
                     "exposure_status": "development_exposed", "model_selection_exposure": "smoke_only", "median_spearman": float(test.spearman.median()),
                     "per_protein_spearman": per_protein_string(test, value_col="spearman"), "random_seeds": "42",
                     "protein_conditionality_pass": bool(pc[pc.shuffle_type == "protein"].prediction_correlation.max() < .995) if not pc.empty else False,
                     "target_conditionality_pass": bool(pc[pc.shuffle_type == "target"].prediction_correlation.max() < .995) if not pc.empty else "not assessed",
                     "evaluation_role": "diagnostic smoke", "comparability_status": "one-fold/not primary", "notes": "Attention entropy near-uniform; not confirmatory"})
    # Dense supervision pilot.
    dense = pd.read_csv(ROOT / "results/v0_5_dense/smoke_results.csv")
    dsum = pd.read_csv(ROOT / "results/v0_5_dense/supervision_delta_summary.csv")
    for (model, protocol), g in dense.groupby(["model", "protocol"]):
        test = g[g.partition == "test"]
        rows.append({"stage": "v0.5-dense", "method": model, "model_family": "pairwise head dense-supervision pilot",
                     "dataset": "GSE237017 designed DBP", "training_proteins": "fold 3 training proteins", "held_out_proteins": "DBP48",
                     "assay": "uPBM processed PBM E-score", "dna_kmer_length": 7, "split_protocol": "one-fold protein_cluster_loco_fold_3",
                     "exposure_status": "development_exposed", "model_selection_exposure": "pilot_fixed", "median_spearman": float(test.spearman.median()),
                     "per_protein_spearman": per_protein_string(test, value_col="spearman"), "random_seeds": "42",
                     "protein_conditionality_pass": "not assessed", "target_conditionality_pass": "not assessed", "evaluation_role": "dense supervision diagnostic",
                     "comparability_status": "one-fold/not primary", "notes": "S512/D4096/D16384; H2 not supported; no scaling continuation"})
    return pd.DataFrame(rows)


def hard_case_analysis():
    base = pd.read_csv(ROOT / "results/v0_5/hard_cases/all_current_models_fail.csv")
    subset_files = {
        "sequence_only_rescued": "results/v0_5/hard_cases/m1_fail_m3_success.csv",
        "conditional_m1c_rescued": "results/v0_5/hard_cases/m1c_fail_m3_success.csv",
        "conditional_m2_rescued": "results/v0_5/hard_cases/m2_fail_m3_success.csv",
        "joint_conditional_rescued": "results/v0_5/hard_cases/joint_controls_fail_m3_success.csv",
        "all_model_failure": "results/v0_5/hard_cases/all_current_models_fail.csv",
    }
    rows = []
    for category, path in subset_files.items():
        d = pd.read_csv(ROOT / path)
        rows.append({"category": category, "count": len(d), "fraction_of_reference_1515": len(d) / 1515,
                     "fraction_of_denominator": len(d) / 1515, "denominator_n": 1515, "reference_population": "1515 preregistered hard cases",
                     "n_proteins": d.dbp_id.nunique(), "proteins": "|".join(sorted(d.dbp_id.unique())),
                     "experimental_score_median": d.experimental_E_score.median(), "target_hamming_median": d.target_hamming.median(),
                     "target_edit_median": d.target_edit.median(), "target_kmer_overlap_median": d.target_kmer_overlap.median(),
                     "gc_fraction_median": d.candidate.map(lambda s: (s.count("G") + s.count("C")) / len(s)).median(),
                     "experimental_percentile_available": False, "interpretation": "pre-registered subset; ranking discrepancy, not binding failure"})
    # Aggregate model sensitivity over frozen inference-only diagnostics.
    for name, file, cond in [("protein_sensitive_cases", "results/v0_6_representation/protein_shuffle.csv", "shuffled_protein"),
                             ("target_sensitive_cases", "results/v0_6_representation/target_shuffle.csv", "shuffled_target")]:
        d = pd.read_csv(ROOT / file)
        corr_col = "prediction_correlation_to_original" if "prediction_correlation_to_original" in d else "prediction_correlation"
        d = d[d.condition == cond] if "condition" in d else d
        rows.append({"category": name, "count": int((d[corr_col] < .995).sum()),
                     "fraction_of_reference_1515": np.nan, "fraction_of_denominator": float((d[corr_col] < .995).mean()), "denominator_n": len(d),
                     "reference_population": "21 v0.6 inference-only shuffle rows", "n_proteins": d.dbp_id.nunique(),
                     "proteins": "|".join(sorted(d.dbp_id.unique())), "experimental_score_median": np.nan,
                     "target_hamming_median": np.nan, "target_edit_median": np.nan, "target_kmer_overlap_median": np.nan,
                     "gc_fraction_median": np.nan, "experimental_percentile_available": False,
                     "interpretation": f"inference-only shuffle sensitivity threshold <0.995; {len(d)} eligible rows"})
    # Descriptive strata over the frozen full candidate-level hard-case replay.
    # Thresholds are fixed before aggregation and are not model-selection gates.
    wide = pd.read_parquet(ROOT / "results/v0_5/phase4_wide_predictions.parquet")
    hard_keys = base[["dbp_id", "candidate"]]
    wide = wide.merge(hard_keys.rename(columns={"candidate": "canonical_7mer"}), on=["dbp_id", "canonical_7mer"], how="inner")
    for category, mask, interpretation in [
        ("motif_near_target", wide["target_kmer_overlap"] >= 0.50, "fixed descriptive stratum: target k-mer overlap >= 0.50"),
        ("highly_dissimilar_off_target", wide["target_hamming"] >= 6 / 7, "fixed descriptive stratum: target Hamming distance >= 6/7"),
        ("high_experimental_low_predicted", (wide["experimental_percentile"] >= 0.90) & (wide["M3_percentile"] <= 0.50), "experimental top decile with M3 at/below median within protein"),
        ("low_experimental_high_predicted", (wide["experimental_percentile"] <= 0.10) & (wide["M3_percentile"] >= 0.90), "experimental bottom decile with M3 at/above 90th percentile"),
    ]:
        d = wide.loc[mask]
        rows.append({"category": category, "count": len(d), "fraction_of_reference_1515": len(d) / 1515,
                     "fraction_of_denominator": len(d) / len(base), "denominator_n": len(base), "reference_population": "871 all-model-failure subset",
                     "n_proteins": d.dbp_id.nunique(), "proteins": "|".join(sorted(d.dbp_id.unique())),
                     "experimental_score_median": d.experimental_score.median(), "target_hamming_median": d.target_hamming.median(),
                     "target_edit_median": d.target_edit.median(), "target_kmer_overlap_median": d.target_kmer_overlap.median(),
                     "gc_fraction_median": d.canonical_7mer.map(lambda s: (s.count("G") + s.count("C")) / len(s)).astype(float).median(),
                     "experimental_percentile_available": True, "interpretation": interpretation})
    out = pd.DataFrame(rows)
    write(out, "hard_case_taxonomy.csv")
    doc = """# Hard Case Analysis\n\nThe reference is the pre-registered v0.3.1 sequence-vs-experiment disagreement set (1,515 candidates). No thresholds were selected after inspecting v0.5/v0.6 predictions. These are ranking discrepancies, not biological binding-failure labels. The `reference_population`, `denominator_n`, `fraction_of_reference_1515` and `fraction_of_denominator` columns distinguish registered 1,515-case subsets, the 871 all-model-failure subset, and inference-only shuffle rows.\n\n""" + md_table(out[["category", "count", "fraction_of_reference_1515", "fraction_of_denominator", "denominator_n", "n_proteins", "proteins", "experimental_score_median", "target_hamming_median", "target_kmer_overlap_median"]]) + """\n\n## Findings\n\nThe strongest reproducible failure class is all-model failure: 871/1,515 (57.5%) candidates remain unresolved by M0-M3. Sequence-only rescue and conditional-model rescue subsets exist, but no stable joint conditional advantage was established (the pre-registered M3 joint-control subset contains 99 candidates). The motif/distance/score strata are calculated within the 871 all-model-failure subset; their denominator is explicitly recorded and they are descriptive only. Existing rows show broad representation across all seven DBPs rather than one protein-specific failure. GC, distance, overlap and score summaries are descriptive only; no small-sample causal interpretation is made. Protein/target sensitivity rows are model-level shuffle diagnostics, not candidate-level biological labels.\n"""
    DOC.mkdir(parents=True, exist_ok=True)
    (DOC / "HARD_CASE_ANALYSIS.md").write_text(doc, encoding="utf-8")


def protein_simulation():
    rng = np.random.default_rng(20260914)
    n_values = [2, 3, 4, 5, 10, 20, 40, 80]
    rows = []
    for rep in range(8):
        zdim, ddim, n_dna, n_test = 8, 12, 256, 20
        W = rng.normal(size=(zdim, ddim))
        dna = rng.normal(size=(n_dna, ddim))
        test_z = rng.normal(size=(n_test, zdim))
        for n_train in n_values:
            train_z = rng.normal(size=(n_train, zdim))
            X = np.concatenate([np.einsum("z,nd->nzd", z, dna).reshape(n_dna, -1) for z in train_z])
            y = np.concatenate([z @ W @ dna.T + rng.normal(0, .5, n_dna) for z in train_z])
            beta = np.linalg.solve(X.T @ X + 1e-2 * np.eye(zdim * ddim), X.T @ y)
            rhos, corrs = [], []
            for i, z in enumerate(test_z):
                pred = np.einsum("nzd,zd->n", np.einsum("z,nd->nzd", z, dna), beta.reshape(zdim, ddim))
                truth = z @ W @ dna.T
                rhos.append(safe_spearman(pred, truth))
                z_shuf = test_z[(i + 1) % n_test]
                shuffled = np.einsum("nzd,zd->n", np.einsum("z,nd->nzd", z_shuf, dna), beta.reshape(zdim, ddim))
                corrs.append(safe_spearman(pred, shuffled))
            rows.append({"simulation": "synthetic_low_rank_protein_specific_landscapes", "replicate": rep,
                         "n_training_proteins": n_train, "n_dna_per_protein": n_dna, "n_test_proteins": n_test,
                         "model_capacity_parameters": zdim * ddim, "held_out_spearman_median": np.nanmedian(rhos),
                         "protein_shuffle_correlation_median": np.nanmedian(corrs), "protein_shuffle_effect_size": 1 - np.nanmedian(corrs),
                         "simulation_only": True, "interpretation": "DNA observations held high; only independent protein units vary"})
    out = pd.DataFrame(rows)
    write(out, "protein_sample_size_simulation.csv")


def natural_designed_shift():
    natural = pd.read_parquet(ROOT / "data/processed/v0_4_1/natural_pbm_benchmark_v0_4_1.parquet")
    designed = pd.read_parquet(ROOT / "data/processed/v0_3_1/designed_dbp_upbm_rc_class_v0_3_1.parquet")
    emb = pd.read_parquet(ROOT / "data/interim/v0_4_2/frozen_plm_embeddings_esm2_t12_35M_UR50D.parquet")
    emb_cols = [c for c in emb.columns if c.startswith("emb_")]
    dseq = designed[["protein_id", "protein_sequence"]].drop_duplicates()
    nseq = natural[["protein_id", "protein_sequence"]].drop_duplicates()
    rows = []
    def aa_comp(seq):
        return np.array([seq.count(a) / len(seq) for a in "ACDEFGHIKLMNPQRSTVWY"])
    for label, df, seqcol, scorecol, dna_col in [("natural", natural, "protein_sequence", "experimental_score", "dna_sequence"), ("designed", designed, "protein_sequence", "experimental_escore_consensus", "canonical_7mer")]:
        g = df.groupby("protein_id")
        for pid, x in g:
            seq = x[seqcol].iloc[0]
            scores = x[scorecol].to_numpy(float)
            top = x.loc[x[scorecol] >= x[scorecol].quantile(.99), dna_col].astype(str).tolist()
            ent = []
            for pos in range(len(top[0])):
                counts = np.bincount(["ACGT".index(s[pos]) for s in top], minlength=4) + 1e-9
                p = counts / counts.sum(); ent.append(float(-(p * np.log2(p)).sum()))
            rows.extend([
                {"metric": "protein_sequence_length", "dataset": label, "protein_id": pid, "value": len(seq), "unit": "aa", "n_units": len(x), "comparability": "descriptive"},
                {"metric": "score_dynamic_range", "dataset": label, "protein_id": pid, "value": float(scores.max() - scores.min()), "unit": "assay_score", "n_units": len(x), "comparability": "assay-specific; not cross-normalized"},
                {"metric": "score_95_5_spread", "dataset": label, "protein_id": pid, "value": float(np.quantile(scores,.95)-np.quantile(scores,.05)), "unit": "assay_score", "n_units": len(x), "comparability": "assay-specific; not cross-normalized"},
                {"metric": "top1pct_position_entropy_mean", "dataset": label, "protein_id": pid, "value": float(np.mean(ent)), "unit": "bits", "n_units": len(top), "comparability": "within-dataset descriptive"},
                {"metric": "dna_kmer_length", "dataset": label, "protein_id": pid, "value": len(str(x[dna_col].iloc[0])), "unit": "nt", "n_units": len(x), "comparability": "protocol difference"},
            ])
    # Sequence composition distribution and cross-domain nearest-neighbor audit.
    nseqs = nseq.set_index("protein_id").protein_sequence.to_dict(); dseqs = dseq.set_index("protein_id").protein_sequence.to_dict()
    for a, seqs in [("natural", nseqs), ("designed", dseqs)]:
        comp = np.vstack([aa_comp(s) for s in seqs.values()])
        for j, aa in enumerate("ACDEFGHIKLMNPQRSTVWY"):
            rows.append({"metric": f"aa_fraction_{aa}", "dataset": a, "protein_id": "__summary__", "value": float(comp[:, j].mean()), "unit": "fraction", "n_units": len(seqs), "comparability": "descriptive"})
    nat_ids = list(nseqs); des_ids = list(dseqs)
    for did in des_ids:
        sims = [SequenceMatcher(None, dseqs[did], nseqs[nid]).ratio() for nid in nat_ids]
        rows.append({"metric": "nearest_natural_sequence_similarity", "dataset": "designed", "protein_id": did, "value": max(sims), "unit": "SequenceMatcher_ratio", "n_units": len(nat_ids), "comparability": "rough sequence similarity; not alignment identity"})
    # Embedding cross-domain distances, if available.
    e = emb.set_index("protein_id")[emb_cols].astype(float)
    for did in des_ids:
        if did not in e.index: continue
        dv = e.loc[did].to_numpy(); vals = []
        for nid in nat_ids:
            if nid not in e.index: continue
            nv = e.loc[nid].to_numpy(); vals.append((1 - np.dot(dv,nv)/(np.linalg.norm(dv)*np.linalg.norm(nv)), np.linalg.norm(dv-nv), nid))
        if vals:
            vals.sort()
            rows.append({"metric": "nearest_natural_embedding_cosine_distance", "dataset": "designed", "protein_id": did, "value": vals[0][0], "unit": "1-cosine", "n_units": len(vals), "comparability": "frozen ESM-2 embedding diagnostic"})
            rows.append({"metric": "nearest_natural_embedding_euclidean_distance", "dataset": "designed", "protein_id": did, "value": vals[0][1], "unit": "L2", "n_units": len(vals), "comparability": "frozen ESM-2 embedding diagnostic"})
    out = pd.DataFrame(rows)
    write(out, "natural_designed_shift.csv")
    agg = out.groupby(["metric", "dataset"], as_index=False).agg(value_median=("value", "median"), value_mean=("value", "mean"), n_rows=("value", "size"))
    doc = """# Natural–Designed Domain Shift\n\nThis is a data-layer comparison, not a new model or transfer run. Natural PBM uses UniPROBE contiguous 8-mers (57 proteins); designed data uses GSE237017 uPBM 7-mers (7 proteins). The local audit reports no verified exact assay constructs for the natural proteins, so protein sequence comparisons are reference-sequence diagnostics.\n\n""" + md_table(agg) + """\n\n## Interpretation\n\nThe two datasets differ in k-mer length, assay source, score processing and construct provenance. Frozen ESM-2 embeddings and rough sequence-similarity nearest neighbors are provided to test whether designed proteins occupy the natural reference distribution, but no PCA/UMAP distance is treated as a formal domain-shift test. Score dynamic range and entropy are kept assay-specific and are not pooled across assays. The correct conclusion is a plausible, measurable domain/construct/assay shift hypothesis, not proof that domain shift alone causes the generalization gap.\n"""
    DOC.mkdir(parents=True, exist_ok=True); (DOC / "NATURAL_DESIGNED_DOMAIN_SHIFT.md").write_text(doc, encoding="utf-8")


def external_and_sota_metadata():
    EXT.mkdir(parents=True, exist_ok=True); META.mkdir(parents=True, exist_ok=True)
    candidates = pd.DataFrame([
        {"study": "GSE237017", "year": 2025, "protein_design_count": 7, "dna_targets": "7-mer uPBM", "assay_type": "uPBM", "full_specificity_landscape": True, "raw_data_available": True, "protein_sequence_available": True, "designed_construct_available": True, "potential_compatibility": "already exposed developmental cohort", "exposure_status": "DEVELOPMENT_EXPOSED; NOT CONFIRMATORY", "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237017"},
        {"study": "PRJNA1014465 (GSE237017-associated sequencing)", "year": 2025, "protein_design_count": "not verified", "dna_targets": "yeast-display sorting / mammalian transcription assays", "assay_type": "sequencing-based functional assay", "full_specificity_landscape": False, "raw_data_available": "SRA metadata indicates deposited reads", "protein_sequence_available": "to verify", "designed_construct_available": "to verify", "potential_compatibility": "orthogonal activity data; not a matched quantitative 7-mer specificity landscape", "exposure_status": "METADATA_ONLY_UNEXPOSED; NOT CONFIRMATORY", "source_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1014465"},
        {"study": "UniPROBE compendium", "year": 2012, "protein_design_count": 57, "dna_targets": "contiguous 8-mer", "assay_type": "PBM", "full_specificity_landscape": True, "raw_data_available": True, "protein_sequence_available": True, "designed_construct_available": False, "potential_compatibility": "natural training/control only; 8-mer and construct mismatch", "exposure_status": "EXPOSED_HISTORICAL_CONTEXT", "source_url": "https://thebrain.bwh.harvard.edu/uniprobe/"},
        {"study": "Prospective new designed-DBP cohort", "year": "TBD", "protein_design_count": "TBD", "dna_targets": "matched 7-mer uPBM preferred", "assay_type": "matched uPBM", "full_specificity_landscape": "to verify", "raw_data_available": "to verify", "protein_sequence_available": "to verify", "designed_construct_available": "to verify", "potential_compatibility": "highest; preregister before labels", "exposure_status": "UNEXPOSED_CONFIRMATORY_CANDIDATE", "source_url": ""},
        {"study": "New public engineered DBP specificity study (candidate search pending)", "year": "TBD", "protein_design_count": "TBD", "dna_targets": "to verify", "assay_type": "to verify", "full_specificity_landscape": "to verify", "raw_data_available": "to verify", "protein_sequence_available": "to verify", "designed_construct_available": "to verify", "potential_compatibility": "metadata-only candidate; no labels read", "exposure_status": "UNEXPOSED_CONFIRMATORY_CANDIDATE", "source_url": ""},
    ])
    candidates.to_csv(EXT / "candidate_datasets.csv", index=False)
    sota = pd.DataFrame([
        {"method": "sequence k-mer / motif similarity", "capability": "sequence-only ranking", "input_requirements": "DNA target/motif", "output_type": "ranking score", "applicable_designed_dbp": True, "full_7mer_landscape": True, "compute": "low", "leakage_risk": "target/motif provenance", "status": "existing comparator; no new run"},
        {"method": "PWM / motif similarity", "capability": "position-specific motif scoring", "input_requirements": "PWM or known motif", "output_type": "PWM score", "applicable_designed_dbp": True, "full_7mer_landscape": True, "compute": "low", "leakage_risk": "motif derived from assay labels", "status": "manifest only"},
        {"method": "SimpleProteinConditionalBaseline", "capability": "protein composition + DNA", "input_requirements": "protein sequence, DNA", "output_type": "ranking score", "applicable_designed_dbp": True, "full_7mer_landscape": True, "compute": "low", "leakage_risk": "natural/design protocol mismatch", "status": "historical unmatched result"},
        {"method": "DeepPBS", "capability": "structure-conditioned PWM ensemble", "input_requirements": "protein-DNA structure / processed structure", "output_type": "PWM/log-probability proxy", "applicable_designed_dbp": "partial", "full_7mer_landscape": True, "compute": "high/Linux", "leakage_risk": "structure overlap/homology", "status": "2/7 diagnostic only"},
        {"method": "NA-MPNN", "capability": "structure-aware sequence design/specificity diagnostic", "input_requirements": "protein-DNA structure", "output_type": "model score", "applicable_designed_dbp": "partial", "full_7mer_landscape": "not established", "compute": "high", "leakage_risk": "PDB overlap and construct mismatch", "status": "not a complete 7-protein baseline"},
        {"method": "structure-conditioned protein-DNA neural models", "capability": "learned contact/sequence compatibility", "input_requirements": "structure or reliable complex prediction", "output_type": "compatibility score", "applicable_designed_dbp": "hypothesis", "full_7mer_landscape": "requires adaptation", "compute": "high", "leakage_risk": "training-set structural homologs", "status": "do not run in v0.7"},
        {"method": "AlphaFold3-derived contact metrics", "capability": "predicted interaction/contact diagnostic", "input_requirements": "protein + DNA complex prediction", "output_type": "contacts/confidence", "applicable_designed_dbp": "hypothesis", "full_7mer_landscape": False, "compute": "very high", "leakage_risk": "model/database overlap; stochasticity", "status": "not a specificity predictor; manifest only"},
        {"method": "Rosetta energetic approaches", "capability": "physics/energy ranking", "input_requirements": "structure and DNA model", "output_type": "energy proxy", "applicable_designed_dbp": "partial", "full_7mer_landscape": "expensive", "compute": "high", "leakage_risk": "structure/template selection", "status": "manifest only"},
    ])
    sota.to_csv(META / "sota_method_manifest.csv", index=False)


def figures_and_docs(master):
    # A tidy, figure-ready long table with no plotted images.
    rows = []
    for _, r in master.iterrows():
        rows.append({"figure": "Figure 1", "panel": "benchmark_design", "stage": r.stage, "method": r.method, "dataset": r.dataset, "split_protocol": r.split_protocol, "exposure_status": r.exposure_status, "median_spearman": r.median_spearman, "conditionality_status": f"P={r.protein_conditionality_pass};T={r.target_conditionality_pass}"})
        if r.stage in ["v0.5", "v0.6", "v0.5-local", "v0.5-dense"]:
            rows.append({"figure": "Figure 2", "panel": "performance_vs_conditionality", "stage": r.stage, "method": r.method, "dataset": r.dataset, "split_protocol": r.split_protocol, "exposure_status": r.exposure_status, "median_spearman": r.median_spearman, "conditionality_status": f"P={r.protein_conditionality_pass};T={r.target_conditionality_pass}"})
    pd.DataFrame(rows).to_csv(OUT / "figure_data_interface.csv", index=False)
    pd.DataFrame([
        {"figure": "Figure 3", "panel": "mechanism_audit", "source": "v0.5 local/dense + v0.6 representation", "status": "available", "data_files": "figure3_mechanism_data.csv", "interpretation": "representation repair changes performance but does not remove real-data conditioning collapse"},
        {"figure": "Figure 4", "panel": "hard_case_and_domain_shift", "source": "v0.7 hard_case_taxonomy + natural_designed_shift", "status": "available", "data_files": "figure4_hard_case_shift_data.csv", "interpretation": "failure taxonomy and shift diagnostics; descriptive, not causal"},
        {"figure": "Figure 5", "panel": "external_validation", "source": "new independent cohort only", "status": "BLOCKED", "data_files": "none", "interpretation": "GSE237017 cannot serve as independent validation"},
    ]).to_csv(OUT / "figure_data_interface_status.csv", index=False)
    # Tidy interfaces for future plotting; values remain directly traceable to
    # frozen artifacts and are not presentation-ready claims.
    mechanism = []
    ra = pd.read_csv(ROOT / "results/v0_6_representation/representation_audit.csv")
    for _, r in ra.iterrows():
        for metric in [c for c in ra.columns if c != "audit"]:
            mechanism.append({"figure": "Figure 3", "panel": "representation_repair", "experiment": r["audit"], "metric": metric, "value": r[metric], "unit": "count_or_boolean", "source": "results/v0_6_representation/representation_audit.csv", "interpretation": "RC-class representation audit"})
    for _, r in pd.read_csv(ROOT / "results/v0_6_representation/go_no_go.csv").iterrows():
        mechanism.append({"figure": "Figure 3", "panel": "conditionality_gates", "experiment": "v0.6_replay", "metric": r["name"], "value": r["observed"], "unit": "gate_observation", "source": "results/v0_6_representation/go_no_go.csv", "interpretation": "pre-registered gate result"})
    for _, r in pd.read_csv(ROOT / "results/v0_5_dense/supervision_delta_summary.csv").iterrows():
        mechanism.append({"figure": "Figure 3", "panel": "dense_supervision", "experiment": r["model"], "metric": "test_spearman_delta", "value": r["D4096_minus_S512"], "unit": "Spearman", "source": "results/v0_5_dense/supervision_delta_summary.csv", "interpretation": "H2 pilot; no scaling claim"})
    for _, r in pd.read_csv(ROOT / "results/v0_5_local/attention_diagnostics.csv").iterrows():
        mechanism.append({"figure": "Figure 3", "panel": "local_attention", "experiment": r["model"], "metric": "attention_entropy_normalized", "value": r["attention_entropy_normalized"], "unit": "normalized_entropy", "source": "results/v0_5_local/attention_diagnostics.csv", "interpretation": "near-uniform attention diagnostic"})
    pd.DataFrame(mechanism).to_csv(OUT / "figure3_mechanism_data.csv", index=False)
    taxonomy = pd.read_csv(OUT / "hard_case_taxonomy.csv")
    shift = pd.read_csv(OUT / "natural_designed_shift.csv")
    f4 = taxonomy.assign(figure="Figure 4", panel="hard_case_taxonomy", source="hard_case_taxonomy.csv", entity=taxonomy.category, metric="count", value=taxonomy["count"], unit="candidates", interpretation="descriptive hard-case stratum")[["figure", "panel", "source", "entity", "metric", "value", "unit", "interpretation"]]
    s4 = shift.groupby(["metric", "dataset"], as_index=False)["value"].median()
    s4 = s4.assign(figure="Figure 4", panel="natural_designed_shift", source="natural_designed_shift.csv", entity=s4["dataset"], unit="dataset_metric", interpretation="descriptive data-layer comparison")
    s4 = s4[["figure", "panel", "source", "entity", "metric", "value", "unit", "interpretation"]]
    pd.concat([f4, s4], ignore_index=True).to_csv(OUT / "figure4_hard_case_shift_data.csv", index=False)
    protocol = """# Benchmark Protocol\n\n## Scope\n\nv0.7 is a frozen evidence synthesis and failure analysis. It reads existing v0.3-v0.6 artifacts and does not train new architectures, search hyperparameters, increase pair counts, or inspect quantitative labels from an external confirmatory cohort.\n\n## Primary axis\n\nThe primary performance unit is within-protein Spearman ranking over canonical reverse-complement DNA units. The v0.5/v0.6 matched benchmark uses the four-fold protein-cluster LOCO split, seeds 17/29/43, the frozen pairwise objective and seven designed proteins.\n\n## Independent axes\n\nPerformance reports median held-out Spearman, per-protein values and seed stability. Conditionality separately reports protein/target shuffle correlation and effect size (1 - correlation). A model with high ranking but shuffle correlation above 0.995 is classified high/low performance plus collapsed conditioning, not successful specificity learning.\n\n## Exposure\n\nAll seven GSE237017 designed DBPs are development-exposed. They are valid for debugging and benchmark diagnosis only. No result in this directory is independent external validation. The master table preserves protocol, exposure and comparability fields; unmatched natural-to-designed SimplePC values remain context, not bridge evidence.\n\n## Reproducibility\n\nHistorical v0.5 files are protected by their frozen manifest. v0.7 output is additive under `results/v0_7_benchmark/`, with machine-readable source paths and status fields.\n"""
    DOC.mkdir(parents=True, exist_ok=True); (DOC / "BENCHMARK_PROTOCOL.md").write_text(protocol, encoding="utf-8")
    assessment = """# V0.7 Scientific Assessment\n\n## Four-state interpretation\n\nThe v0.5/v0.6 conditional models occupy **low-to-moderate performance + collapsed conditioning** on the exposed designed replay: v0.6 all-seven medians are M0 0.2829, M1 0.2016, M1c 0.1693, M2 0.1342 and M3 0.1415, while protein-shuffle median correlation is 1.0000 and target-shuffle median correlation is 0.9969. The sequence-only comparator is 0.2321. This prevents interpreting raw Spearman as proof of protein-specific specificity.\n\n## Evidence and limits\n\nThere is sufficient evidence for a **benchmark-observed generalization gap under this strict unseen-protein protocol**: the gap is replicated across v0.5, v0.6, local/residue and dense-supervision diagnostics, and 871/1,515 hard cases remain unresolved by all current models. It is not yet sufficient to claim a universal biological impossibility. Seven independent designed proteins make a small-N artifact plausible; the simulation explicitly demonstrates that large DNA N cannot substitute for protein-level diversity, but it is not an empirical threshold.\n\n## Formal claims now supported\n\n1. Under the current GSE237017 design, strict unseen-protein ranking is modest and conditional heads do not show measurable protein sensitivity.\n2. v0.5 RC averaging caused severe cross-class representation collisions; v0.6 repair improves developmental ranking but does not remove collapse.\n3. Dense pair scaling and residue/local attention pilots did not establish a conditional advantage.\n4. Natural PBM and designed uPBM differ in k-mer length, assay, score semantics and construct provenance; existing SimplePC transfer numbers are unmatched context.\n\n## Hypotheses only\n\nNatural-to-designed domain shift, assay/construct mismatch, insufficient protein diversity and inadequate sequence representations may each contribute. Their causal contributions cannot be separated with seven exposed proteins and unmatched natural constructs.\n\n## Route decision\n\n**Route B: benchmark / generalization-gap paper**, with a data-expansion prerequisite. Route A is not supported because no current conditional method passes real-data conditionality. Route C remains a requirement for a stronger causal claim: obtain a new independent designed-DBP cohort and preregister the frozen evaluation before reading labels. Do not present GSE237017 as Figure 5 external validation.\n\n## Critical missing experiment\n\nThe single highest-value experiment is a genuinely independent, assay-matched designed-DBP cohort with substantially more than seven protein units, full 7-mer specificity landscapes, construct sequences and a preregistered frozen benchmark. This simultaneously tests whether the observed gap survives exposure control and whether protein-level sample size, rather than DNA-level pair count, is limiting.\n\n## Strongest claim\n\n**Under strict unseen-protein evaluation, existing sequence-based approaches show a reproducible benchmark-level generalization gap on the currently available de novo DBP cohort: representation repair improves ranking but does not restore protein- or target-dependent prediction, and the evidence is not yet sufficient for a new predictive-method claim or independent validation claim.**\n"""
    (DOC / "V0_7_SCIENTIFIC_ASSESSMENT.md").write_text(assessment, encoding="utf-8")


def build_performance_conditionality(master):
    """Join performance rows to only the shuffle diagnostics that exist."""
    rows = []
    seq_baseline = float(master.loc[master.method.eq("sequence_kmer3"), "median_spearman"].iloc[0])
    v06 = pd.read_csv(ROOT / "results/v0_6_representation/conditionality_diagnostics.csv").set_index("model")
    v05p = pd.read_csv(ROOT / "results/v0_5/shuffled_protein_diagnostic.csv")
    v05t = pd.read_csv(ROOT / "results/v0_5/shuffled_target_diagnostic.csv")
    v06p = pd.read_csv(ROOT / "results/v0_6_representation/protein_shuffle.csv")
    v06t = pd.read_csv(ROOT / "results/v0_6_representation/target_shuffle.csv")
    local_sh = pd.read_csv(ROOT / "results/v0_5_local/shuffle_diagnostics.csv")
    for _, r in master.iterrows():
        stage, method = r.stage, r.method
        p_corr = t_corr = p_eff = t_eff = np.nan
        p_status = t_status = "not_assessed"
        if stage == "v0.5" and method == "M3":
            p_corr = v05p.loc[v05p.condition.eq("shuffled_protein"), "prediction_correlation_to_original"].median()
            t_corr = v05t.loc[v05t.condition.eq("shuffled_target"), "prediction_correlation_to_original"].median()
            p_eff, t_eff = 1 - p_corr, 1 - t_corr
        elif stage == "v0.6" and method == "M3":
            p_corr = v06p.prediction_correlation.median(); t_corr = v06t.prediction_correlation.median()
            p_eff, t_eff = 1 - p_corr, 1 - t_corr
        elif stage == "v0.5-local":
            p = local_sh[(local_sh.model == method) & (local_sh.shuffle_type == "protein")]
            t = local_sh[(local_sh.model == method) & (local_sh.shuffle_type == "target")]
            if not p.empty: p_corr, p_eff = p.prediction_correlation.median(), 1 - p.prediction_correlation.median()
            if not t.empty: t_corr, t_eff = t.prediction_correlation.median(), 1 - t.prediction_correlation.median()
        if pd.notna(p_corr): p_status = "collapsed" if p_corr > .995 else "real_conditioning"
        if pd.notna(t_corr): t_status = "collapsed" if t_corr > .995 else "real_conditioning"
        perf = "high" if float(r.median_spearman) >= seq_baseline else "low"
        cond = "collapsed" if p_status == "collapsed" or t_status == "collapsed" else ("real_conditioning" if "real_conditioning" in {p_status, t_status} else "not_assessed")
        rows.append({"stage": stage, "method": method, "dataset": r.dataset, "split_protocol": r.split_protocol,
                     "exposure_status": r.exposure_status, "median_heldout_spearman": r.median_spearman,
                     "per_protein_spearman": r.per_protein_spearman, "random_seeds": r.random_seeds,
                     "protein_shuffle_correlation": p_corr, "target_shuffle_correlation": t_corr,
                     "protein_conditioning_effect_size": p_eff, "target_conditioning_effect_size": t_eff,
                     "protein_conditionality_status": p_status, "target_conditionality_status": t_status,
                     "performance_axis": perf, "conditionality_axis": cond,
                     "four_state": f"{perf}_performance+{cond}_conditioning", "high_performance_reference": f"sequence-only median {seq_baseline:.6f}",
                     "collapse_threshold": .995, "notes": "NA means diagnostic not run; no inference from missing shuffle"})
    write(pd.DataFrame(rows), "performance_conditionality.csv")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    master = build_master(); master.to_csv(OUT / "benchmark_master.csv", index=False)
    hard_case_analysis(); protein_simulation(); natural_designed_shift(); external_and_sota_metadata(); figures_and_docs(master)
    build_performance_conditionality(master)
    # A concise provenance manifest for automated checks.
    pd.DataFrame([{"artifact": p.name, "source": "v0.7 generated", "quantitative_external_labels_read": False} for p in OUT.glob("*.csv")]).to_csv(OUT / "v0_7_artifact_manifest.csv", index=False)


if __name__ == "__main__":
    main()
