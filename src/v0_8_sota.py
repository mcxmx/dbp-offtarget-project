"""v0.8 frozen SOTA audit and external-validation firewall.

This module only reads frozen v0.3-v0.7 artifacts.  It does not train, tune,
or inspect quantitative labels from a new external cohort.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DBPS = ["DBP1", "DBP3", "DBP35", "DBP48", "DBP5", "DBP6", "DBP9"]
COLS = ["protein_id", "canonical_7mer"]


def spearman(a: pd.Series, b: pd.Series) -> float:
    x = pd.Series(a).rank(method="average").to_numpy(float)
    y = pd.Series(b).rank(method="average").to_numpy(float)
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def ensure_dirs() -> dict[str, Path]:
    dirs = {
        "results": ROOT / "results/v0_8_sota",
        "sota_meta": ROOT / "metadata/v0_8_sota",
        "external_meta": ROOT / "metadata/v0_8_external",
        "sota_docs": ROOT / "docs/v0_8_sota",
        "external_docs": ROOT / "docs/v0_8_external",
        "paper": ROOT / "docs/paper",
        "tests": ROOT / "tests/v0_8",
    }
    for p in dirs.values():
        p.mkdir(parents=True, exist_ok=True)
    return dirs


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def standardized_predictions() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    exp = pd.read_parquet(ROOT / "data/processed/v0_3_1/designed_dbp_upbm_rc_class_v0_3_1.parquet")
    exp = exp[COLS + ["experimental_escore_consensus"]].rename(
        columns={"experimental_escore_consensus": "experimental_score"}
    )
    deep = pd.read_parquet(ROOT / "results/v0_4_2/tables/deeppbs_predictions_completed_v0_4_2.parquet")
    deep = deep.rename(columns={"deeppbs_score": "prediction_score"})
    deep["method"] = "DeepPBS"
    deep["scoring_protocol_frozen"] = "max contiguous 7-mer window/orientation; sum log(P+1e-9)"
    na = pd.read_parquet(ROOT / "results/v0_4/tables/nampnn_predictions.parquet")
    na = na.rename(columns={"prediction_score": "prediction_score"})
    na["method"] = "NA-MPNN"
    na["scoring_protocol_frozen"] = "official PPM best-window log probability; same canonical RC convention"
    for df in (deep, na):
        df["exposure_status"] = "DEVELOPMENT_EXPOSED_GSE237017"
        df["evaluation_role"] = "partial structure diagnostic; not independent validation"
    deep_out = deep[["protein_id", "canonical_7mer", "prediction_score", "method", "structure_id",
                     "structure_source", "model_version", "prediction_type", "scoring_protocol_frozen",
                     "exposure_status", "evaluation_role"]].copy()
    na_out = na[["protein_id", "canonical_7mer", "prediction_score", "method", "structure_id",
                 "model_version", "prediction_type", "scoring_protocol_frozen", "exposure_status",
                 "evaluation_role"]].copy()
    return exp, deep_out, na_out


def method_metrics(exp: pd.DataFrame, pred: pd.DataFrame, method: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    joined = exp.merge(pred[COLS + ["prediction_score"]], on=COLS, how="left")
    rows = []
    for protein in DBPS:
        sub = joined[joined.protein_id == protein].dropna(subset=["prediction_score", "experimental_score"])
        rows.append({
            "method": method,
            "protein_id": protein,
            "status": "evaluated" if len(sub) else "not_evaluable_missing_structure_prediction",
            "n_candidates": int(len(sub)),
            "spearman": spearman(sub.experimental_score, sub.prediction_score) if len(sub) else np.nan,
            "coverage_fraction": float(len(sub) / 8192),
            "score_variance": float(sub.prediction_score.var()) if len(sub) else np.nan,
            "exposure_status": "DEVELOPMENT_EXPOSED_GSE237017",
        })
    per = pd.DataFrame(rows)
    eligible = per.loc[per.status == "evaluated", "spearman"].dropna()
    summary = pd.DataFrame([{
        "method": method,
        "metric": "spearman",
        "median_spearman": float(eligible.median()) if len(eligible) else np.nan,
        "mean_spearman": float(eligible.mean()) if len(eligible) else np.nan,
        "n_proteins_evaluated": int(len(eligible)),
        "coverage": f"{int(len(eligible))}/7",
        "protein_unit_n": int(len(eligible)),
        "statistical_note": "protein is the independent unit; no candidate-level significance",
    }])
    return per, summary


def hard_case_metrics(pred: pd.DataFrame, method: str) -> dict[str, object]:
    hard = pd.read_parquet(ROOT / "results/v0_4_2/tables/baseline_failure_cases_deeppbs_completed_v0_4_2.parquet")
    p = pred[COLS + ["prediction_score"]].copy()
    p["predicted_percentile"] = p.groupby("protein_id")["prediction_score"].rank(pct=True)
    h = hard.merge(p, on=COLS, how="left")
    eligible = h.dropna(subset=["predicted_percentile"])
    resolved = eligible[eligible.predicted_percentile >= 0.90]
    return {
        "method": method,
        "hard_case_reference_n": int(len(hard)),
        "hard_case_evaluable_n": int(len(eligible)),
        "hard_case_resolved_n": int(len(resolved)),
        "hard_case_resolution_rate": float(len(resolved) / len(eligible)) if len(eligible) else np.nan,
        "all_model_failure_reference_n": 871,
        "note": "resolution means predicted percentile >= 0.90; exposed diagnostic only",
    }


def parse_per_protein(path: Path, model_cols: list[str], stage: str) -> list[dict[str, object]]:
    df = pd.read_csv(path)
    out = []
    for model in model_cols:
        vals = df[["dbp_id", model]].dropna()
        vals = vals.set_index("dbp_id")[model]
        out.append({
            "stage": stage,
            "method": model,
            "dataset": "GSE237017 designed DBP",
            "protocol_class": "protein_cluster_loco",
            "training_proteins": "fold-specific subset of 7",
            "held_out_proteins": "7 fold-held-out rows",
            "assay": "uPBM processed E-score",
            "dna_kmer_length": 7,
            "exposure_status": "DEVELOPMENT_EXPOSED_ALL7",
            "model_selection_exposure": "frozen historical result; no v0.8 selection",
            "median_spearman": float(vals.median()),
            "per_protein_spearman": "|".join(f"{k}:{v:.6f}" for k, v in vals.items()),
            "protein_unit_n": int(len(vals)),
            "random_seeds": "17|29|43",
            "protein_conditioning_evidence": "collapsed for M3; not assessed for M0-M2",
            "target_conditioning_evidence": "collapsed for M3; not assessed for M0-M2",
            "comparability_status": "matched historical benchmark",
        })
    return out


def create_master(per_struct: pd.DataFrame, struct_summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    # Historical sequence-only baseline.
    seq = pd.read_csv(ROOT / "results/v0_3_1/tables/designed_dbp_sequence_baseline_rc_aware.csv")
    s = seq[seq.metric.eq("kmer3_jaccard_to_paper_motif_rc_aware")] if "metric" in seq else seq
    if len(s):
        vals = s.set_index("protein_id")["spearman"] if "spearman" in s else s.set_index("protein_id")["metric_value"]
        rows.append({"method": "sequence_kmer3", "model_family": "sequence-only", "protein_conditioned": False,
                     "target_conditioned": False, "requires_complex_structure": False, "full_landscape": True,
                     "coverage": "7/7", "median_spearman": float(vals.median()),
                     "per_protein_spearman": "|".join(f"{k}:{v:.6f}" for k, v in vals.items()),
                     "protein_unit_n": int(len(vals)), "hard_case_recovery": "historical v0.7 taxonomy",
                     "protocol_class": "v0.3.1 historical sequence baseline", "exposure_status": "DEVELOPMENT_EXPOSED_ALL7",
                     "comparability_status": "sequence proxy; matched candidate universe", "method_version": "v0.3.1 frozen",
                     "code_repository": "local repository", "checkpoint_version": "none", "required_input": "target motif",
                     "output_type": "k-mer Jaccard ranking", "runtime": "historical", "failed_cases": "see v0.7 taxonomy"})
    # v0.5 frozen and v0.6 corrected replay rows.
    rows.extend({"method": r["method"], "model_family": "sequence conditional pairwise head", "protein_conditioned": r["method"] in {"M1", "M1c", "M3"},
                 "target_conditioned": r["method"] in {"M2", "M3"}, "requires_complex_structure": False, "full_landscape": True,
                 "coverage": "7/7", "median_spearman": r["median_spearman"], "per_protein_spearman": r["per_protein_spearman"],
                 "protein_unit_n": 7, "hard_case_recovery": "see v0.7 hard_case_taxonomy.csv", "protocol_class": r["protocol_class"],
                 "exposure_status": r["exposure_status"], "comparability_status": r["comparability_status"], "method_version": r.get("stage", "v0.6 corrected"),
                 "code_repository": "local repository", "checkpoint_version": "frozen historical checkpoint", "required_input": "protein/candidate/target per M0-M3",
                 "output_type": "ranking score", "runtime": "historical", "failed_cases": "see v0.7 taxonomy"}
                for r in parse_per_protein(ROOT / "results/v0_5/primary_per_protein_results.csv", ["M0", "M1", "M1c", "M2", "M3"], "v0.5")
                + [{"method": r["method"], "model_family": "corrected sequence conditional pairwise head", "protein_conditioned": r["method"] in {"M1", "M1c", "M3"},
                    "target_conditioned": r["method"] in {"M2", "M3"}, "requires_complex_structure": False, "full_landscape": True,
                    "coverage": "7/7", "median_spearman": r["median_spearman"], "per_protein_spearman": r["per_protein_spearman"],
                    "protein_unit_n": 7, "hard_case_recovery": "see v0.7 hard_case_taxonomy.csv", "protocol_class": "v0.6 corrected RC + ordered target replay",
                    "exposure_status": "DEVELOPMENT_EXPOSED_ALL7", "comparability_status": "matched replay; representation repair only", "method_version": "v0.6 corrected",
                    "code_repository": "local repository", "checkpoint_version": "seeds 17|29|43", "required_input": "protein/target/candidate per M0-M3",
                    "output_type": "ranking score", "runtime": "historical replay", "failed_cases": "see v0.7 taxonomy"}
                   for r in parse_per_protein(ROOT / "results/v0_6_representation/v06_per_protein_spearman.csv", ["M0", "M1", "M1c", "M2", "M3"], "v0.6")])
    for _, r in struct_summary.iterrows():
        p = per_struct[(per_struct.method == r.method) & (per_struct.status == "evaluated")]
        rows.append({"method": r.method, "model_family": "structure-conditioned specificity proxy", "protein_conditioned": True,
                     "target_conditioned": False, "requires_complex_structure": True, "full_landscape": True,
                     "coverage": r.coverage, "median_spearman": r.median_spearman,
                     "per_protein_spearman": "|".join(f"{x.protein_id}:{x.spearman:.6f}" for _, x in p.iterrows()),
                     "protein_unit_n": r.n_proteins_evaluated, "hard_case_recovery": "see hard_case_taxonomy.csv",
                     "protocol_class": "v0.4 frozen structure diagnostic; PWM/PPM proxy", "exposure_status": "DEVELOPMENT_EXPOSED_ALL7",
                     "comparability_status": "partial coverage; not a 7-protein generalization estimate", "method_version": "frozen v0.4.2",
                     "code_repository": "official public repository", "checkpoint_version": "official bundled ensemble / s_70114.pt", "required_input": "protein-DNA complex structure",
                     "output_type": "PWM/PPM log-probability proxy", "runtime": "historical Linux/CPU", "failed_cases": "5/7 no compatible structure"})
    rows += [
        {"method": "SimpleProteinConditionalBaseline", "model_family": "composition + ridge", "protein_conditioned": True, "target_conditioned": False,
         "requires_complex_structure": False, "full_landscape": True, "coverage": "7/7", "median_spearman": 0.36161081294792224,
         "per_protein_spearman": "historical artifact; see v0_4_1", "protein_unit_n": 7, "hard_case_recovery": "not re-run",
         "protocol_class": "natural-trained unmatched transfer", "exposure_status": "DEVELOPMENT_EXPOSED_ALL7", "comparability_status": "context only; unmatched natural/design protocol", "method_version": "v0.4.1", "code_repository": "local repository", "checkpoint_version": "ridge fit", "required_input": "protein sequence + DNA", "output_type": "ranking score", "runtime": "historical", "failed_cases": "not re-run"},
        {"method": "FrozenPLMProteinConditionalBaseline", "model_family": "frozen ESM-2 + ridge", "protein_conditioned": True, "target_conditioned": False,
         "requires_complex_structure": False, "full_landscape": True, "coverage": "7/7", "median_spearman": 0.1530196201433334,
         "per_protein_spearman": "historical artifact; see v0_4_2", "protein_unit_n": 7, "hard_case_recovery": "historical v0.4.2",
         "protocol_class": "natural-trained unmatched transfer", "exposure_status": "DEVELOPMENT_EXPOSED_ALL7", "comparability_status": "context only; unmatched natural/design protocol", "method_version": "v0.4.2", "code_repository": "local repository", "checkpoint_version": "frozen ESM-2", "required_input": "protein sequence + DNA", "output_type": "ranking score", "runtime": "historical", "failed_cases": "not re-run"},
        {"method": "experimental_replicate_reference", "model_family": "experiment vs experiment", "protein_conditioned": False, "target_conditioned": False,
         "requires_complex_structure": False, "full_landscape": True, "coverage": "7/7", "median_spearman": 0.5914386964424022,
         "per_protein_spearman": "assay replicate median; see experimental_reference.csv", "protein_unit_n": 7, "hard_case_recovery": "not applicable",
         "protocol_class": "assay reproducibility context", "exposure_status": "DEVELOPMENT_EXPOSED_ALL7", "comparability_status": "not a model and not a mathematical ceiling", "method_version": "v0.3.1", "code_repository": "local repository", "checkpoint_version": "none", "required_input": "replicate assay measurements", "output_type": "Spearman reference", "runtime": "historical", "failed_cases": "not applicable"},
        {"method": "AF3_contact_probability", "model_family": "predicted complex contact diagnostic", "protein_conditioned": True, "target_conditioned": True,
         "requires_complex_structure": True, "full_landscape": False, "coverage": "0/7; feasibility only", "median_spearman": np.nan,
         "per_protein_spearman": "not run", "protein_unit_n": 0, "hard_case_recovery": "not run",
         "protocol_class": "pre-registered subset structural diagnostic", "exposure_status": "DEVELOPMENT_EXPOSED_ALL7", "comparability_status": "not a full-landscape specificity predictor", "method_version": "not run", "code_repository": "official AF3", "checkpoint_version": "not available", "required_input": "protein + DNA complex", "output_type": "contact probability", "runtime": "not run", "failed_cases": "57,344 full inferences infeasible"},
        {"method": "ContactSeek_contact_probability", "model_family": "AF3 confidence/contact metric", "protein_conditioned": True, "target_conditioned": True,
         "requires_complex_structure": True, "full_landscape": False, "coverage": "0/7; exposed example only", "median_spearman": np.nan,
         "per_protein_spearman": "not run", "protein_unit_n": 0, "hard_case_recovery": "not run",
         "protocol_class": "subset contact diagnostic", "exposure_status": "DEVELOPMENT_EXPOSED_ALL7", "comparability_status": "not affinity and not independent validation", "method_version": "public workflow; not run", "code_repository": "ContactSeek", "checkpoint_version": "AF3-dependent", "required_input": "AF3 confidence outputs", "output_type": "contact probability", "runtime": "not run", "failed_cases": "exposed example only"},
    ]
    # These rows preserve prior negative diagnostics without implying a new run.
    rows.extend([
        {"method": "local_residue_attention", "model_family": "local/residue attention pilot", "protein_conditioned": True, "target_conditioned": True, "requires_complex_structure": False, "full_landscape": False, "coverage": "pilot", "median_spearman": np.nan, "per_protein_spearman": "see results/v0_5_local", "protein_unit_n": 2, "hard_case_recovery": "not comparable", "protocol_class": "v0.5 diagnostic pilot", "exposure_status": "DEVELOPMENT_EXPOSED", "comparability_status": "attention near-uniform; no rescue", "method_version": "v0.5", "code_repository": "local repository", "checkpoint_version": "pilot", "required_input": "protein residues + DNA", "output_type": "ranking score", "runtime": "historical pilot", "failed_cases": "conditioning collapsed"},
        {"method": "dense_supervision", "model_family": "dense pair supervision pilot", "protein_conditioned": True, "target_conditioned": True, "requires_complex_structure": False, "full_landscape": False, "coverage": "pilot", "median_spearman": np.nan, "per_protein_spearman": "see results/v0_5_dense", "protein_unit_n": 1, "hard_case_recovery": "not comparable", "protocol_class": "v0.5 diagnostic pilot", "exposure_status": "DEVELOPMENT_EXPOSED", "comparability_status": "pair scaling did not rescue collapse", "method_version": "v0.5", "code_repository": "local repository", "checkpoint_version": "pilot", "required_input": "protein/target/candidate pairs", "output_type": "ranking score", "runtime": "historical pilot", "failed_cases": "conditioning collapsed"},
    ])
    return pd.DataFrame(rows)


def write_docs(dirs: dict[str, Path], deep_summary: pd.DataFrame, na_summary: pd.DataFrame) -> None:
    (dirs["sota_docs"] / "SOTA_BENCHMARK_PROTOCOL.md").write_text("""# V0.8 SOTA Benchmark Protocol\n\nThis is a frozen audit of already-run structure diagnostics plus historical sequence baselines. No model is fine-tuned on GSE237017, no hyperparameter is selected from the exposed cohort, and no new sequence architecture is introduced.\n\n## Primary unit\n\nThe independent statistical unit is protein (N=7 for the designed cohort). Candidate-level scores (8,192 canonical RC classes per protein) are used only to calculate within-protein Spearman and ranking metrics. Structure methods are reported with their actual coverage (N=2 for DeepPBS and NA-MPNN), never as a seven-protein estimate.\n\n## Uniform scoring\n\nDeepPBS PWM and NA-MPNN PPM outputs are converted to a fixed canonical 7-mer log-probability ranking. No PBM E-score numerical equivalence is claimed. Register/window and strand conventions are inherited from the frozen v0.4 artifacts. Missing structures are missing, not imputed.\n\n## Conditionality\n\nProtein/target shuffle diagnostics are reported only where the model defines a counterfactual. Structure scores are not assigned meaningless shuffle values because each structure is protein-specific. A high Spearman with no conditionality evidence is classified as collapsed, not successful specificity learning.\n\n## Exposure\n\nAll seven GSE237017 DBPs are development-exposed. These results are benchmark/development diagnostics and cannot be called external validation.\n""", encoding="utf-8")
    (dirs["sota_docs"] / "AF3_FEASIBILITY_AUDIT.md").write_text("""# AF3 Feasibility Audit\n\nA full 7-protein x 8,192-candidate complex landscape would require 57,344 separate protein-DNA complex predictions before contact extraction, with additional sensitivity to stochastic inference and DNA register. AF3 contact probability is a structural/contact diagnostic, not a calibrated affinity or PBM E-score predictor. ContactSeek's public workflow extracts contact probabilities from AF3 confidence JSON/embedding outputs and demonstrates an exposed Glasscock-designed example; it does not provide a frozen, independent seven-protein full landscape for this project.\n\nNo defensible cheap approximation was accepted: replacing per-candidate complex inference with a sequence-only surrogate would change the estimand and violate the structure-method audit. `af3_candidate_manifest.csv` therefore records a preregistered subset definition only (on-target, Hamming-1/2, the frozen hard-case set, and fixed-seed stratified background); `af3_contact_metrics.csv` is an explicit not-run schema. No AF3 result was used for method selection.\n\nStatus: `FEASIBILITY_ONLY`, not a performance benchmark and not external validation.\n""", encoding="utf-8")
    (dirs["sota_docs"] / "V0_8_SOTA_RESULTS.md").write_text(f"""# V0.8 SOTA Results\n\n## Structure diagnostics\n\nDeepPBS: median Spearman {deep_summary.iloc[0].median_spearman:.6f} across {int(deep_summary.iloc[0].n_proteins_evaluated)}/7 proteins (DBP35 and DBP48 only). NA-MPNN is evaluated on the same two structure-covered proteins; its exact metrics are in `na_mpnn_metrics.csv`. Neither result is a seven-protein transfer estimate.\n\nThe experimental replicate reference is approximately 0.5914 Spearman (experiment-vs-experiment), a context reference rather than a mathematical ceiling. Sequence-only and v0.6 conditional rows remain exposed developmental comparators.\n\n## Interpretation\n\nThe available structure rows do not establish that structure methods close the gap because coverage is only N=2 and both structures are development-exposed. AF3/ContactSeek feasibility does not create a full-landscape score. The strongest defensible claim remains limited to sequence-based approaches on the current cohort.\n\n## Scientific judgment\n\n1. The present evidence supports a reproducible sequence-model benchmark gap, not a universal computational impossibility.\n2. DeepPBS and NA-MPNN do not currently provide evidence of a rescue; their partial medians are below the replicate reference and are underpowered.\n3. The 871/1515 common failures cannot be declared structure-method failures because only a subset is structure-evaluable.\n4. The paper title should say **sequence-based specificity prediction** unless a future, adequately covered structure benchmark changes the result.\n5. No independent confirmatory cohort is currently qualified without label access and preregistration; the project is **BLOCKED ON INDEPENDENT CONFIRMATORY DATA**.\n6. Route B remains the appropriate paper route, with data expansion as the key blocker.\n""", encoding="utf-8")
    (dirs["external_docs"] / "EXTERNAL_COHORT_SELECTION_PROTOCOL.md").write_text("""# External Cohort Selection Protocol\n\nThis protocol is written before quantitative label access. A confirmatory cohort must contain de novo designed DNA-binding proteins, construct/protein sequences, intended targets, quantitative specificity measurements, multiple independently designed proteins (preferably >7), and an assay compatible with the current canonical 7-mer ranking task. Complete 7-mer landscapes or a sufficiently dense randomized library are preferred.\n\nMetadata inspection may establish study existence, assay type, sequence availability and raw-data availability. It may not inspect per-sequence scores, performance labels, easy/hard protein annotations, or any result used for model selection. GSE237017 and ContactSeek's Glasscock example are exposed and excluded from confirmation.\n\nA cohort is not \"qualified\" until all inclusion fields are verified without labels, the frozen validation plan is committed, and the manifest is hashed. If no cohort meets the criteria, report `NO SUITABLE PUBLIC CONFIRMATORY COHORT IDENTIFIED` and request data from authors.\n""", encoding="utf-8")
    (dirs["external_docs"] / "PREREGISTERED_VALIDATION_PLAN.md").write_text("""# Preregistered Validation Plan (Template; Labels Locked)\n\n**Status:** FROZEN TEMPLATE, not executed. No external quantitative labels have been read.\n\n1. **Primary hypothesis:** frozen methods have measurable transferable DNA ranking on an independent designed-DBP cohort under strict unseen-protein evaluation.\n2. **Methods:** sequence k-mer, frozen v0.6 M0-M3, SimplePC/FrozenPLM context only, and any structure method with verified compatible inputs; no fine-tuning.\n3. **Checkpoint/code:** hashes are recorded in `FROZEN_VALIDATION_MANIFEST.yaml`.\n4. **Scoring:** canonical 7-mer within-protein ranking; fixed PWM/PPM log-probability conversion for structure methods.\n5. **Register/strand:** predeclared biological register; reverse-complement classes are scored invariantly; ambiguous register means exclusion, never best-result alignment.\n6. **Primary metric:** median per-protein Spearman.\n7. **Secondary metrics:** per-protein Spearman, top-k/NDCG where pre-specified, hard-case recovery, and conditionality where meaningful.\n8. **Exclusions:** missing construct/target, incompatible assay, ambiguous register, or incomplete scores according to predeclared thresholds.\n9. **Missing data:** retain all eligible proteins in denominator; report coverage and failures.\n10. **Failure handling:** preserve failed outputs and logs; no silent retries with changed settings.\n11. **No post-hoc selection:** no model, checkpoint, register, metric, or protein removal after label access.\n12. **Statistics:** protein is the independent unit; paired exact/permutation summaries and effect sizes; no candidate-level significance claims.\n13. **Confirmatory vs exploratory:** the primary frozen pipeline is confirmatory; every later analysis is labeled POST-HOC / EXPLORATORY.\n\nExternal validation remains blocked until a cohort satisfies the selection protocol.\n""", encoding="utf-8")
    (dirs["paper"] / "PAPER_OUTLINE_V0_8.md").write_text("""# Paper Outline V0.8\n\n## Introduction\nDefine transfer of DNA-binding specificity to unseen de novo proteins and distinguish performance from conditionality.\n\n## Result 1: Strict unseen-protein benchmark\n- Claim: current sequence-based models show modest transfer and a reproducible benchmark gap.\n- Supporting files: `results/v0_7_benchmark/benchmark_master.csv`, `results/v0_8_sota/sota_benchmark_master.csv`.\n- Figure: Figure 1, cohort/split schematic.\n- Limitation: seven exposed designed proteins.\n- Reviewer objection: protein-level N is small and methods are heterogeneous.\n\n## Result 2: Representation repair\n- Claim: RC collision repair improves ranking but does not restore conditionality.\n- Supporting files: `docs/v0_6_representation/`, `results/v0_6_representation/`.\n- Figure: Figure 2, performance vs conditionality.\n- Limitation: replay remains development-exposed.\n- Reviewer objection: repair is necessary but not sufficient.\n\n## Result 3: Failed rescue attempts\n- Claim: dense supervision and local/residue attention did not rescue conditioning collapse.\n- Supporting files: `results/v0_5_dense/`, `results/v0_5_local/`.\n- Figure: Figure 3.\n- Limitation: pilots are not universal architecture proof.\n- Reviewer objection: other structure-aware methods were not yet fully covered.\n\n## Result 4: Structure-method audit\n- Claim: DeepPBS/NA-MPNN are legal partial diagnostics but current coverage is N=2; AF3 contact metrics are feasibility-only.\n- Supporting files: `docs/v0_8_sota/`, `results/v0_8_sota/`.\n- Figure: coverage-aware benchmark matrix.\n- Limitation: cannot infer universal structure-method failure.\n- Reviewer objection: request broader structure coverage.\n\n## Result 5: Hard cases and domain shift\n- Claim: 871/1515 common sequence-model failures and measurable natural/designed assay/construct differences motivate, but do not prove, a domain-shift/data-regime explanation.\n- Supporting files: v0.7 hard taxonomy and shift CSV.\n- Figure: Figure 4.\n- Limitation: causal contributions are unidentified.\n- Reviewer objection: hard cases are label/noise dependent.\n\n## Result 6: Independent validation\n- Status: blocked; do not include a Figure 5 until a preregistered cohort is available.\n\n## Discussion / Limitations / Methods\nEmphasize protein-level N=7, development exposure, no universal computational claim, and the highest-value next experiment: an independent assay-matched cohort with >7 designed proteins and dense specificity landscapes.\n""", encoding="utf-8")


def create_external_metadata(dirs: dict[str, Path]) -> None:
    rows = [
        {"study": "GSE237017 / Glasscock et al.", "year": 2025, "protein_design_count": 7, "dna_targets": "7-mer uPBM", "assay_type": "uPBM", "full_specificity_landscape": True, "raw_data_available": True, "protein_sequence_available": True, "designed_construct_available": True, "potential_compatibility": "assay compatible but already exposed", "exposure_status": "DEVELOPMENT_EXPOSED_NOT_CONFIRMATORY", "label_access": False, "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237017"},
        {"study": "RFdiffusion3 DNA-binding protein study", "year": 2026, "protein_design_count": "15 target campaigns; exact qualifying cohort to verify", "dna_targets": "target-specific plus randomized-library assays reported in metadata", "assay_type": "variant competition / randomized library screening", "full_specificity_landscape": "to verify", "raw_data_available": "to verify", "protein_sequence_available": "to verify", "designed_construct_available": "to verify", "potential_compatibility": "promising metadata-only candidate; assay/register compatibility unverified", "exposure_status": "UNEXPOSED_CONFIRMATORY_CANDIDATE", "label_access": False, "source_url": "https://www.biorxiv.org/content/10.64898/2026.04.27.720408v1"},
        {"study": "ContactSeek Glasscock example", "year": 2025, "protein_design_count": "1 exposed example", "dna_targets": "on-target and mismatch examples", "assay_type": "AF3 contact probability", "full_specificity_landscape": False, "raw_data_available": True, "protein_sequence_available": True, "designed_construct_available": True, "potential_compatibility": "contact diagnostic only; exposed", "exposure_status": "DEVELOPMENT_EXPOSED_NOT_CONFIRMATORY", "label_access": False, "source_url": "https://github.com/menghaowei/ContactSeek"},
        {"study": "PRJNA1014465 / GSE237017-associated", "year": 2025, "protein_design_count": "not verified", "dna_targets": "functional assay sequencing", "assay_type": "sequencing-based activity", "full_specificity_landscape": False, "raw_data_available": True, "protein_sequence_available": "to verify", "designed_construct_available": "to verify", "potential_compatibility": "orthogonal activity; not matched 7-mer landscape", "exposure_status": "METADATA_ONLY_UNEXPOSED_NOT_CONFIRMATORY", "label_access": False, "source_url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1014465"},
    ]
    write_csv(pd.DataFrame(rows), dirs["external_meta"] / "candidate_datasets.csv")
    write_csv(pd.DataFrame([{"target": "RFdiffusion3 authors", "request": "exact designed construct sequences, intended DNA targets, randomized-library data, processed quantitative tables, assay metadata", "status": "REQUEST_REQUIRED_BEFORE_CONFIRMATION"}, {"target": "new engineered DBP study authors", "request": "same matched construct and specificity package", "status": "REQUEST_REQUIRED_BEFORE_CONFIRMATION"}]), dirs["external_meta"] / "data_request_targets.csv")
    (dirs["external_docs"] / "DATA_REQUEST_EMAIL.md").write_text("""# Data Request Email\n\nSubject: Request for designed DNA-binding protein specificity data\n\nDear authors,\n\nWe are preparing a preregistered benchmark of DNA-binding specificity transfer to unseen de novo designed proteins. Could you share, if available, the exact protein/construct sequences, intended target DNA sequences, randomized-library or specificity-assay data, processed quantitative tables, and assay/normalization metadata for the designed constructs?\n\nThe requested data would be used only after a frozen analysis plan is registered. We will not use it to select models or tune parameters, and we will report failures and missing data.\n\nBest regards,\nThe DBP specificity benchmark team\n""", encoding="utf-8")


def create_manifests(dirs: dict[str, Path]) -> None:
    rows = [
        {"method": "DeepPBS", "version": "official commit 8bfb211dd67f02877841f6f33aa493ddf7daedf", "code_repository": "https://github.com/timkartar/DeepPBS", "checkpoint_version": "official bundled ensemble", "required_input": "protein-DNA complex or processed structure", "structure_source": "DBP035.pdb and 8TAC helix-only", "target_dna_requirements": "DNA geometry; fixed 7-mer window scoring", "output_type": "PWM/log-probability proxy", "runtime": "completed historical Linux run", "coverage": "2/7", "failed_cases": "5 proteins missing reliable compatible structures", "manual_intervention": "documented v0.4.2 preprocessing", "exposure_leakage": "GSE237017 development-exposed; no exact overlap found in checked manifests; homolog unknown", "status": "FROZEN_PARTIAL_DIAGNOSTIC"},
        {"method": "NA-MPNN", "version": "specificity checkpoint s_70114.pt", "code_repository": "https://github.com/baker-laboratory/NA-MPNN", "checkpoint_version": "s_70114.pt", "required_input": "protein-DNA complex PDB/CIF and DNA positions", "structure_source": "DBP035.pdb and 8TAC helix-only", "target_dna_requirements": "DNA mask/window; fixed canonical RC 7-mer scoring", "output_type": "PPM/log-probability proxy", "runtime": "completed historical Linux/CPU-compatible run", "coverage": "2/7", "failed_cases": "not a complete seven-protein run", "manual_intervention": "official specificity preprocessing", "exposure_leakage": "GSE237017 development-exposed; DBP48 overlap caveat retained", "status": "FROZEN_PARTIAL_DIAGNOSTIC"},
        {"method": "AlphaFold3 contact probability", "version": "not run", "code_repository": "https://github.com/google-deepmind/alphafold3", "checkpoint_version": "not available locally", "required_input": "protein + DNA complex inference per candidate", "structure_source": "predicted", "target_dna_requirements": "candidate-specific DNA", "output_type": "contact/confidence diagnostic", "runtime": "not run; feasibility only", "coverage": "0/7", "failed_cases": "cost/availability and estimand mismatch", "manual_intervention": "none", "exposure_leakage": "subset would use exposed GSE237017 only", "status": "FEASIBILITY_ONLY"},
        {"method": "ContactSeek", "version": "public workflow; not run", "code_repository": "https://github.com/menghaowei/ContactSeek", "checkpoint_version": "AF3 confidence JSON/embedding dependent", "required_input": "AF3 predicted complex outputs", "structure_source": "predicted", "target_dna_requirements": "on/off-target candidate subset", "output_type": "contact probability / delta contact", "runtime": "not run", "coverage": "0/7 full landscape", "failed_cases": "not an affinity predictor; exposed example", "manual_intervention": "none", "exposure_leakage": "Glasscock example is development-exposed", "status": "SUBSET_DIAGNOSTIC_NOT_RUN"},
    ]
    write_csv(pd.DataFrame(rows), dirs["sota_meta"] / "method_execution_manifest.csv")
    write_csv(pd.DataFrame([{"protein_id": p, "candidate_set": "on_target; Hamming-1; Hamming-2; frozen 1515 hard cases; 50 fixed-seed background", "selection_seed": 202608, "defined_before_af3": True, "status": "NOT_RUN_FEASIBILITY_ONLY"} for p in DBPS]), dirs["results"] / "af3_candidate_manifest.csv")
    write_csv(pd.DataFrame(columns=["protein_id", "canonical_7mer", "contact_probability", "delta_contact_probability", "af3_model_version", "status"]), dirs["results"] / "af3_contact_metrics.csv")


def main() -> None:
    dirs = ensure_dirs()
    exp, deep, na = standardized_predictions()
    deep.to_csv(dirs["results"] / "deeppbs_predictions.csv", index=False)
    na.to_csv(dirs["results"] / "na_mpnn_predictions.csv", index=False)
    deep_per, deep_sum = method_metrics(exp, deep, "DeepPBS")
    na_per, na_sum = method_metrics(exp, na, "NA-MPNN")
    deep_per.to_csv(dirs["results"] / "deeppbs_metrics.csv", index=False)
    na_per.to_csv(dirs["results"] / "na_mpnn_metrics.csv", index=False)
    hard_rows = [hard_case_metrics(deep, "DeepPBS"), hard_case_metrics(na, "NA-MPNN")]
    hard = pd.read_csv(ROOT / "results/v0_7_benchmark/hard_case_taxonomy.csv")
    hard = pd.concat([hard, pd.DataFrame([
        {"category": "DeepPBS_structure_evaluable", "count": hard_rows[0]["hard_case_evaluable_n"], "fraction_of_reference_1515": hard_rows[0]["hard_case_evaluable_n"] / 1515, "fraction_of_denominator": 1.0, "denominator_n": hard_rows[0]["hard_case_evaluable_n"], "reference_population": "1515 hard cases; eligible structure subset", "n_proteins": 2, "proteins": "DBP35|DBP48", "experimental_score_median": np.nan, "target_hamming_median": np.nan, "target_edit_median": np.nan, "target_kmer_overlap_median": np.nan, "gc_fraction_median": np.nan, "experimental_percentile_available": True, "interpretation": "structure coverage, not universal structure-method failure"},
        {"category": "NA-MPNN_structure_evaluable", "count": hard_rows[1]["hard_case_evaluable_n"], "fraction_of_reference_1515": hard_rows[1]["hard_case_evaluable_n"] / 1515, "fraction_of_denominator": 1.0, "denominator_n": hard_rows[1]["hard_case_evaluable_n"], "reference_population": "1515 hard cases; eligible structure subset", "n_proteins": 2, "proteins": "DBP35|DBP48", "experimental_score_median": np.nan, "target_hamming_median": np.nan, "target_edit_median": np.nan, "target_kmer_overlap_median": np.nan, "gc_fraction_median": np.nan, "experimental_percentile_available": True, "interpretation": "structure coverage, not universal structure-method failure"},
    ])], ignore_index=True)
    hard.to_csv(dirs["results"] / "hard_case_taxonomy.csv", index=False)
    # Existing v0.7 simulation and data-layer shift are copied into a new immutable stage.
    pd.read_csv(ROOT / "results/v0_7_benchmark/protein_sample_size_simulation.csv").assign(stage="v0.8_reused_simulation").to_csv(dirs["results"] / "protein_sample_size_simulation.csv", index=False)
    pd.read_csv(ROOT / "results/v0_7_benchmark/natural_designed_shift.csv").assign(stage="v0.8_reused_data_layer_audit").to_csv(dirs["results"] / "natural_designed_shift.csv", index=False)
    write_csv(pd.read_csv(ROOT / "results/v0_7_benchmark/benchmark_master.csv"), dirs["results"] / "benchmark_master_v0_7_replay.csv")
    master = create_master(pd.concat([deep_per.assign(method="DeepPBS"), na_per.assign(method="NA-MPNN")], ignore_index=True), pd.concat([deep_sum, na_sum], ignore_index=True))
    for h in hard_rows:
        mask = master.method.eq(h["method"])
        master.loc[mask, "hard_case_recovery"] = f"{h['hard_case_resolved_n']}/{h['hard_case_evaluable_n']} ({h['hard_case_resolution_rate']:.3f}); exposed diagnostic"
    master.to_csv(dirs["results"] / "sota_benchmark_master.csv", index=False)
    cond = pd.DataFrame([
        {"method": "v0.6 M3", "performance_median_spearman": 0.1415054997760777, "protein_shuffle_correlation": 1.0, "target_shuffle_correlation": 0.9969, "protein_conditioning_effect_size": 0.0, "target_conditioning_effect_size": 0.0031, "performance_axis": "low", "conditionality_axis": "collapsed", "four_state": "low_performance+collapsed_conditioning", "protein_unit_n": 7},
        {"method": "DeepPBS", "performance_median_spearman": float(deep_sum.iloc[0].median_spearman), "protein_shuffle_correlation": np.nan, "target_shuffle_correlation": np.nan, "protein_conditioning_effect_size": np.nan, "target_conditioning_effect_size": np.nan, "performance_axis": "low", "conditionality_axis": "not_applicable_structure_specific", "four_state": "low_performance+not_assessed", "protein_unit_n": 2},
        {"method": "NA-MPNN", "performance_median_spearman": float(na_sum.iloc[0].median_spearman), "protein_shuffle_correlation": np.nan, "target_shuffle_correlation": np.nan, "protein_conditioning_effect_size": np.nan, "target_conditioning_effect_size": np.nan, "performance_axis": "low", "conditionality_axis": "not_applicable_structure_specific", "four_state": "low_performance+not_assessed", "protein_unit_n": 2},
    ])
    cond.to_csv(dirs["results"] / "performance_conditionality.csv", index=False)
    rep = pd.read_csv(ROOT / "results/v0_3_1/tables/experimental_noise_ceiling.csv")
    rep["reference_type"] = "experiment_vs_experiment"
    rep["independent_unit"] = "protein"
    rep["interpretation"] = "empirical assay reproducibility reference; not a strict mathematical ceiling"
    rep.to_csv(dirs["results"] / "experimental_reference.csv", index=False)
    create_manifests(dirs)
    create_external_metadata(dirs)
    write_docs(dirs, deep_sum, na_sum)
    # Hash the frozen v0.8 registration inputs after generation.
    tracked = [ROOT / "src/v0_8_sota.py", dirs["external_docs"] / "EXTERNAL_COHORT_SELECTION_PROTOCOL.md", dirs["external_docs"] / "PREREGISTERED_VALIDATION_PLAN.md"]
    hashes = {str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in tracked}
    for p in [ROOT / "external/nampnn/NA-MPNN/models/specificity_model/s_70114.pt", ROOT / "external/deeppbs/DeepPBS/run/plot_scripts/txts/DeepPBS.txt"]:
        hashes[str(p.relative_to(ROOT)).replace("\\", "/")] = hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else "MISSING"
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False).stdout.strip()
    manifest = {"status": "pre_external_template_not_active", "label_access": False, "git_commit": commit, "git_tag_proposal": "v0.8-pre-external-validation", "date": "2026-09-14", "code_and_checkpoint_sha256": hashes, "dataset_metadata_known_before_label_access": ["candidate_datasets.csv", "EXTERNAL_COHORT_SELECTION_PROTOCOL.md"], "no_label_access_assertion": True}
    import yaml
    (dirs["external_meta"] / "FROZEN_VALIDATION_MANIFEST.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    print(json.dumps({"deepPBS_median": float(deep_sum.iloc[0].median_spearman), "NA_MPNN_median": float(na_sum.iloc[0].median_spearman), "deepPBS_coverage": int(deep_sum.iloc[0].n_proteins_evaluated), "NA_MPNN_coverage": int(na_sum.iloc[0].n_proteins_evaluated)}, indent=2))


if __name__ == "__main__":
    main()
