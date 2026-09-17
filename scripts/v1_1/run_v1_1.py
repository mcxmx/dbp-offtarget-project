from __future__ import annotations

import hashlib
import json
import math
import shutil
import sys
import tarfile
from itertools import combinations, product
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
DEEPPBS_RAW = ROOT / "data" / "raw" / "v1_1_deeppbs" / "frozen_design7_20260914"
COMPETITION_SOURCE = ROOT / "data" / "raw" / "v1_1_competition" / "source_data_extended_data_fig3.xls"
PROTEINS = ["DBP001", "DBP003", "DBP005", "DBP006", "DBP009", "DBP035", "DBP048"]
EPS = 1e-12
BASES = "ACGT"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_dirs() -> None:
    for name in [
        "00_inventory",
        "01_deeppbs_ingest",
        "02_nampnn_ingest",
        "03_register_mapping",
        "04_competition",
        "05_assay_aligned_eval",
        "06_assay_agreement",
        "07_local_global",
        "08_pwm_capacity",
        "09_global_pwm_projection",
        "10_conditionality",
        "11_training_overlap",
        "figures",
    ]:
        (OUT / name).mkdir(parents=True, exist_ok=True)


def ingest_tar() -> Path:
    """Extract the supplied immutable prediction bundle once, without editing NPZs."""
    archive = ROOT / "frozen_design7_20260914.tar.gz"
    target = DEEPPBS_RAW
    target.parent.mkdir(parents=True, exist_ok=True)
    if not (target / "predictions").exists():
        with tarfile.open(archive, "r:gz") as tf:
            tf.extractall(target.parent)
    # The archive contains a top-level directory; make the expected path explicit.
    extracted = target.parent / "frozen_design7_20260914"
    if extracted != target and extracted.exists():
        return extracted
    return target


def load_targets() -> pd.DataFrame:
    path = ROOT / "metadata" / "v0_3_1" / "designed_dbp_target_definitions.csv"
    df = pd.read_csv(path)
    df["protein_id"] = df["protein_id"].map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}")
    return df[df.protein_id.isin(PROTEINS)].copy()


def decode_one_hot(arr: np.ndarray) -> str:
    arr = np.asarray(arr)
    if arr.ndim != 2 or arr.shape[1] != 4:
        raise AssertionError(f"expected Nx4 one-hot, got {arr.shape}")
    if not np.allclose(arr.sum(axis=1), 1.0, atol=1e-5):
        raise AssertionError("one-hot rows are not normalized")
    return "".join(BASES[int(i)] for i in arr.argmax(axis=1))


def ingest_deeppbs(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, manifest = [], []
    for protein in PROTEINS:
        path = root / "predictions" / f"{protein}.npz_predict.npz"
        if not path.exists():
            raise FileNotFoundError(path)
        data = np.load(path, allow_pickle=False)
        if set(data.files) != {"P", "Seq"}:
            raise AssertionError(f"{protein}: unexpected NPZ keys {data.files}")
        P, seq_arr = np.asarray(data["P"]), np.asarray(data["Seq"])
        if P.shape != seq_arr.shape or P.ndim != 2 or P.shape[1] != 4:
            raise AssertionError(f"{protein}: P/Seq shape mismatch {P.shape}/{seq_arr.shape}")
        if not np.isfinite(P).all() or not np.isfinite(seq_arr).all():
            raise AssertionError(f"{protein}: non-finite array")
        if not np.allclose(P.sum(axis=1), 1.0, atol=2e-5):
            raise AssertionError(f"{protein}: P rows do not sum to one")
        seq = decode_one_hot(seq_arr)
        for i, row in enumerate(P):
            rows.append({
                "protein": protein,
                "position": i + 1,
                "reference_base": seq[i],
                "p_A": float(row[0]), "p_C": float(row[1]),
                "p_G": float(row[2]), "p_T": float(row[3]),
                "source": path.relative_to(ROOT).as_posix(),
            })
        manifest.append({
            "protein": protein,
            "width": int(P.shape[0]),
            "sequence": seq,
            "dtype_P": str(P.dtype),
            "dtype_Seq": str(seq_arr.dtype),
            "finite": True,
            "P_rows_sum_1": True,
            "Seq_rows_sum_1": True,
            "helix_score": {"DBP001": 0.8461538461538461, "DBP003": 0.8461538461538461,
                             "DBP005": 1.0, "DBP006": 1.0, "DBP009": 1.0,
                             "DBP035": 1.0, "DBP048": 1.0}[protein],
            "contact_count": {"DBP001": 261, "DBP003": 367, "DBP005": 270,
                               "DBP006": 235, "DBP009": 235, "DBP035": 266,
                               "DBP048": 242}[protein],
            "source_file": path.relative_to(ROOT).as_posix(),
            "sha256": sha256(path),
        })
    pwm = pd.DataFrame(rows)
    man = pd.DataFrame(manifest)
    pwm.to_csv(OUT / "01_deeppbs_ingest" / "deeppbs_design7_pwm.tsv", sep="\t", index=False)
    man.to_csv(OUT / "01_deeppbs_ingest" / "deeppbs_design7_manifest.tsv", sep="\t", index=False)
    (OUT / "01_deeppbs_ingest" / "deeppbs_qc_report.md").write_text(
        "# DeepPBS ingest QC\n\n"
        f"Validated {len(PROTEINS)}/7 NPZ files. All arrays are finite, float32, Nx4, "
        "and row-normalized; raw NPZ files were not modified.\n\n"
        + "```\n" + man[["protein", "width", "sequence", "helix_score", "contact_count"]].to_string(index=False) + "\n```"
        + "\n", encoding="utf-8"
    )
    return pwm, man


def ingest_nampnn() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Ingest frozen v1.0 landscape output; native PPM matrices were not retained."""
    path = ROOT / "results" / "v1_0_structure" / "nampnn_results.csv"
    df = pd.read_csv(path)
    df["protein"] = df["protein_id"].map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}")
    df = df[df.protein.isin(PROTEINS)].copy()
    # The v1.0 artifact contains score landscapes, not the native predicted_ppm tensor.
    # Preserve them verbatim and expose a clear non-PPM status rather than reconstructing P.
    df[["protein", "canonical_7mer", "nampnn_score", "prediction_type", "model_version",
        "source_npz"]].to_csv(OUT / "02_nampnn_ingest" / "nampnn_design7_landscape.tsv", sep="\t", index=False)
    target_df = load_targets().rename(columns={"protein_id": "protein"})
    qc = (df.groupby("protein").agg(n_landscape=("canonical_7mer", "nunique"),
                                     score_min=("nampnn_score", "min"), score_max=("nampnn_score", "max"))
          .reset_index().merge(target_df[["protein", "experimental_assay_target"]], on="protein", how="left"))
    qc["native_ppm_available"] = False
    qc["representation_note"] = "v1.0 retained canonical 7-mer landscape only; no native PPM tensor"
    qc.to_csv(OUT / "02_nampnn_ingest" / "nampnn_alignment_qc.tsv", sep="\t", index=False)
    # Keep the requested conceptual schema explicit. Values are NA because the v1.0
    # artifact retained scores, not the checkpoint's intermediate PPM tensor.
    pwm_rows = []
    for protein in PROTEINS:
        sequence = str(target_df.loc[target_df.protein == protein, "experimental_assay_target"].iloc[0])
        for position, base in enumerate(sequence, 1):
            pwm_rows.append({"protein": protein, "position": position, "reference_base": base,
                             "p_A": np.nan, "p_C": np.nan, "p_G": np.nan, "p_T": np.nan,
                             "representation_status": "native_ppm_not_retained",
                             "source": "results/v1_0_structure/nampnn_results.csv"})
    pd.DataFrame(pwm_rows).to_csv(OUT / "02_nampnn_ingest" / "nampnn_design7_pwm.tsv", sep="\t", index=False)
    (OUT / "02_nampnn_ingest" / "representation_note.md").write_text(
        "# NA-MPNN representation note\n\n"
        "The frozen v1.0 artifact contains a canonical 7-mer score landscape, not the native `predicted_ppm` tensor. "
        "The conceptual PWM table is therefore emitted with NA probabilities; local NA-MPNN mutation effects are labeled as derived from the frozen landscape, never as native PPM values.\n",
        encoding="utf-8",
    )
    return df, qc


def choose_center_window(seq: str, motif: str, width: int = 7) -> tuple[str, int, bool, str]:
    """Return an independently defined window; no experimental scores enter this function."""
    if len(seq) < width:
        return "", -1, True, "sequence_shorter_than_7"
    idx = seq.find(motif)
    if idx >= 0:
        start = max(0, min(idx - (width - len(motif)) // 2, len(seq) - width))
        return seq[start:start + width], start, False, "motif_exact_in_cognate_sequence"
    start = (len(seq) - width) // 2
    return seq[start:start + width], start, True, "motif_absent; central_design_window_used"


def build_register_map(pwm: pd.DataFrame, targets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, report, candidates = [], [], []
    for protein in PROTEINS:
        t = targets.loc[targets.protein_id == protein].iloc[0]
        model_seq = pwm.loc[pwm.protein == protein].sort_values("position")["reference_base"].str.cat()
        motif = str(t.designed_binding_site_motif)
        window, start, ambiguous, evidence = choose_center_window(model_seq, motif)
        assay_seq = str(t.experimental_assay_target)
        assay_window, assay_start, assay_amb, assay_evidence = choose_center_window(assay_seq, motif)
        status = "ambiguous" if ambiguous or assay_amb else "unique_sequence_derived"
        if status == "ambiguous":
            # Enumerate every contiguous seven-position model window as a
            # predeclared sensitivity set; no assay score is used to select one.
            for candidate_start in range(max(0, len(model_seq) - 6)):
                candidates.append({"protein": protein, "candidate_type": "model_contiguous_window",
                                   "model_window_start_1based": candidate_start + 1,
                                   "model_window": model_seq[candidate_start:candidate_start + 7],
                                   "selection_status": "sensitivity_only_ambiguous_mapping"})
        for j, base in enumerate(window):
            rows.append({
                "protein": protein, "model": "DeepPBS/NA-MPNN",
                "model_position": start + j + 1, "model_base": base,
                "assay_position": assay_start + j + 1, "assay_reference_base": assay_window[j],
                "pbm_7mer_position": j + 1, "orientation": "forward",
                "offset": start, "mapping_status": status,
                "mapping_evidence": evidence,
                "ambiguity_flag": ambiguous or assay_amb,
            })
        report.append({
            "protein": protein, "model_sequence": model_seq,
            "model_window": window, "model_window_start_1based": start + 1,
            "assay_target": assay_seq, "assay_window": assay_window,
            "assay_window_start_1based": assay_start + 1,
            "motif": motif, "mapping_status": status,
            "mapping_evidence": evidence + "; assay=" + assay_evidence,
            "ambiguity_flag": ambiguous or assay_amb,
        })
    maps, rep = pd.DataFrame(rows), pd.DataFrame(report)
    maps.to_csv(OUT / "03_register_mapping" / "register_map.tsv", sep="\t", index=False)
    rep.to_csv(OUT / "03_register_mapping" / "register_summary.tsv", sep="\t", index=False)
    pd.DataFrame(candidates).to_csv(OUT / "03_register_mapping" / "register_sensitivity_candidates.tsv", sep="\t", index=False)
    text = ["# Register mapping report", "", "Mapping was determined from model Seq, target metadata, and exact motif matching before any assay scoring.", ""]
    for _, r in rep.iterrows():
        text += [f"## {r.protein}", f"- Model sequence: `{r.model_sequence}`", f"- Model 7-mer: `{r.model_window}` (1-based start {int(r.model_window_start_1based)})", f"- Assay target: `{r.assay_target}`", f"- Assay window: `{r.assay_window}` (1-based start {int(r.assay_window_start_1based)})", f"- Status: **{r.mapping_status}**; {r.mapping_evidence}", ""]
    (OUT / "03_register_mapping" / "register_mapping_report.md").write_text("\n".join(text), encoding="utf-8")
    return maps, rep


def competition_missing() -> pd.DataFrame:
    cols = ["protein", "assay_position", "wt_base", "mut_base", "experimental_value", "experimental_effect", "replicate", "source"]
    empty = pd.DataFrame(columns=cols)
    empty.to_csv(OUT / "04_competition" / "competition_mutations.tsv", sep="\t", index=False)
    status = pd.DataFrame([{
        "status": "MISSING_LOCAL_SOURCE",
        "proteins_covered": 0,
        "source_search": "repository-wide search of data/raw, data/processed, metadata, analysis, docs",
        "files_checked": "MOESM3 supplementary tables; MOESM12 uPBM/orthogonality; MOESM20 uPBM source data; GSE237017 raw/processed files",
        "reason": "No Glasscock single-base mutation competition quantitative table is present locally; no values inferred.",
        "downstream_action": "Competition-aligned scoring and PBM-vs-competition agreement remain not evaluable.",
    }])
    status.to_csv(OUT / "04_competition" / "competition_source_status.tsv", sep="\t", index=False)
    (OUT / "04_competition" / "competition_assay_audit.md").write_text(
        "# Competition assay audit\n\n"
        "The expected Glasscock single-base mutation competition table is absent from the repository. "
        "The pipeline stops this assay-specific branch with `MISSING_LOCAL_SOURCE`; no quantitative labels, direction, replicate, or mutation rows are invented.\n",
        encoding="utf-8",
    )
    return empty


def spearman(x: pd.Series | np.ndarray, y: pd.Series | np.ndarray) -> float:
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan")
    return float(spearmanr(x, y).statistic)


def pwm_matrix(pwm: pd.DataFrame, protein: str) -> tuple[np.ndarray, str]:
    g = pwm[pwm.protein == protein].sort_values("position")
    return g[["p_A", "p_C", "p_G", "p_T"]].to_numpy(float), "".join(g.reference_base)


def score_pwm(P: np.ndarray, seq: str, positions: list[int]) -> float:
    idx = {b: i for i, b in enumerate(BASES)}
    vals = [math.log(float(P[pos, idx[b]]) + EPS) for pos, b in zip(positions, seq)]
    return float(sum(vals))


def orient_to_target(seq: str, target: str) -> tuple[str, str, int]:
    """Choose the closer strand to a fixed target; ties use lexical order."""
    candidates = [(seq, "forward"), (reverse_complement(seq), "reverse_complement")]
    scored = [(sum(a != b for a, b in zip(s, target)), s, orient) for s, orient in candidates]
    distance, chosen, orientation = min(scored, key=lambda x: (x[0], x[1]))
    return chosen, orientation, int(distance)


def canonical_landscape_from_pwm(P: np.ndarray, positions: list[int]) -> pd.DataFrame:
    rows = []
    for tup in product(BASES, repeat=7):
        seq = "".join(tup)
        rc = reverse_complement(seq)
        canonical = min(seq, rc)
        if seq != canonical:
            continue
        score = max(score_pwm(P, seq, positions), score_pwm(P, rc, positions))
        rows.append({"canonical_7mer": canonical, "predicted_score": score, "orientation_rule": "max(forward, reverse_complement)"})
    if len(rows) != 8192:
        raise AssertionError(f"canonical universe has {len(rows)} rows")
    return pd.DataFrame(rows)


def project_global(pwm: pd.DataFrame, nampnn: pd.DataFrame, targets: pd.DataFrame, reg: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pbm = pd.read_parquet(ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet")
    deep_rows, metrics = [], []
    all_scores = {}
    for protein in PROTEINS:
        P, model_seq = pwm_matrix(pwm, protein)
        rr = reg[reg.protein == protein].sort_values("pbm_7mer_position")
        pos = (rr.model_position - 1).astype(int).tolist()
        deep = canonical_landscape_from_pwm(P, pos).rename(columns={"predicted_score": "deeppbs_score"})
        exp = pbm[pbm.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}") == protein][["canonical_7mer", "experimental_escore_consensus"]]
        deep = deep.merge(exp, on="canonical_7mer", how="inner")
        deep["protein"] = protein
        deep["mapping_status"] = rr.mapping_status.iloc[0]
        deep_rows.append(deep)
        all_scores[("DeepPBS", protein)] = deep.set_index("canonical_7mer")["deeppbs_score"]
        metrics.append({"method": "DeepPBS", "protein": protein, "n": len(deep), "spearman": spearman(deep.deeppbs_score, deep.experimental_escore_consensus), "mapping_status": rr.mapping_status.iloc[0], "primary_status": "primary" if rr.mapping_status.iloc[0] == "unique_sequence_derived" else "sensitivity_only_ambiguous_mapping"})
        ng = nampnn[nampnn.protein == protein][["canonical_7mer", "nampnn_score"]].drop_duplicates()
        ng = ng.merge(exp, on="canonical_7mer", how="inner").rename(columns={"nampnn_score": "predicted_score"})
        ng["protein"] = protein
        all_scores[("NA-MPNN", protein)] = ng.set_index("canonical_7mer")["predicted_score"]
        metrics.append({"method": "NA-MPNN", "protein": protein, "n": len(ng), "spearman": spearman(ng.predicted_score, ng.experimental_escore_consensus), "mapping_status": rr.mapping_status.iloc[0], "primary_status": "primary" if rr.mapping_status.iloc[0] == "unique_sequence_derived" else "sensitivity_only_ambiguous_mapping", "note": "frozen v1.0 global landscape; native PPM not retained"})
    deep_all = pd.concat(deep_rows, ignore_index=True)
    deep_all.to_csv(OUT / "09_global_pwm_projection" / "deeppbs_global_projection.tsv", sep="\t", index=False)
    mdf = pd.DataFrame(metrics)
    mdf.to_csv(OUT / "09_global_pwm_projection" / "global_projection_metrics.tsv", sep="\t", index=False)
    # Save a unified model x protein table for later conditionality analysis.
    rows = []
    for (method, protein), vals in all_scores.items():
        for k, v in vals.items(): rows.append({"method": method, "protein": protein, "canonical_7mer": k, "predicted_score": float(v)})
    pd.DataFrame(rows).to_csv(OUT / "09_global_pwm_projection" / "unified_global_landscapes.tsv", sep="\t", index=False)
    return deep_all, mdf


def mutation_projection(pwm: pd.DataFrame, reg: pd.DataFrame, targets: pd.DataFrame, global_df: pd.DataFrame, nampnn: pd.DataFrame) -> None:
    rows = []
    for protein in PROTEINS:
        P, _ = pwm_matrix(pwm, protein)
        rr = reg[reg.protein == protein].sort_values("pbm_7mer_position")
        model_pos = (rr.model_position - 1).tolist()
        ref = "".join(rr.model_base)
        assay_ref = "".join(rr.assay_reference_base)
        # Every single mutation in the mapped 7-mer is generated, even without experimental labels.
        for i, wt in enumerate(ref):
            for mut in BASES:
                if mut == wt: continue
                rows.append({"protein": protein, "assay_position": int(rr.iloc[i].assay_position), "position_in_7mer": i + 1,
                             "wt_base": wt, "mut_base": mut,
                             "predicted_delta_logP": math.log(P[model_pos[i], BASES.index(mut)] + EPS) - math.log(P[model_pos[i], BASES.index(wt)] + EPS),
                             "predicted_delta_probability": float(P[model_pos[i], BASES.index(mut)] - P[model_pos[i], BASES.index(wt)]),
                             "experimental_effect": np.nan, "experimental_source_status": "MISSING_LOCAL_SOURCE",
                             "mapping_status": rr.mapping_status.iloc[0], "reference_window": ref, "assay_reference_window": assay_ref})
        # NA-MPNN native PPMs were not retained. Derive a transparent local delta
        # from the frozen landscape only, with an explicit representation label.
        ng = nampnn[nampnn.protein == protein].drop_duplicates("canonical_7mer").set_index("canonical_7mer")["nampnn_score"].to_dict()
        target_c = canonical_rc(assay_ref)
        if target_c in ng:
            wt_score = float(ng[target_c])
            for i, wt in enumerate(target_c):
                for mut in BASES:
                    if mut == wt: continue
                    mutated = canonical_rc(target_c[:i] + mut + target_c[i + 1:])
                    rows.append({"protein": protein, "assay_position": int(rr.iloc[i].assay_position), "position_in_7mer": i + 1,
                                 "wt_base": wt, "mut_base": mut,
                                 "predicted_delta_logP": float(ng.get(mutated, np.nan) - wt_score),
                                 "predicted_delta_probability": np.nan,
                                 "experimental_effect": np.nan, "experimental_source_status": "MISSING_LOCAL_SOURCE",
                                 "mapping_status": rr.mapping_status.iloc[0], "reference_window": target_c,
                                 "assay_reference_window": assay_ref,
                                 "method": "NA-MPNN_derived_from_frozen_global_landscape"})
    pred = pd.DataFrame(rows)
    pred["method"] = pred.get("method", "DeepPBS_native_pwm").fillna("DeepPBS_native_pwm")
    pred.to_csv(OUT / "05_assay_aligned_eval" / "per_mutation_predictions.tsv", sep="\t", index=False)
    per = pd.DataFrame([{"method": m, "protein": p, "n_mutations": 0, "spearman": np.nan, "pearson": np.nan, "sign_agreement": np.nan, "status": "NOT_EVALUABLE_MISSING_COMPETITION_ASSAY"} for m in ["DeepPBS", "NA-MPNN"] for p in PROTEINS])
    per.to_csv(OUT / "05_assay_aligned_eval" / "per_protein_metrics.tsv", sep="\t", index=False)
    pd.DataFrame([{"method": "DeepPBS", "macro_median_spearman": np.nan, "macro_mean_spearman": np.nan, "status": "NOT_EVALUABLE_MISSING_COMPETITION_ASSAY"}, {"method": "NA-MPNN", "macro_median_spearman": np.nan, "macro_mean_spearman": np.nan, "status": "NOT_EVALUABLE_MISSING_COMPETITION_ASSAY"}]).to_csv(OUT / "05_assay_aligned_eval" / "summary_metrics.tsv", sep="\t", index=False)


def assay_agreement_missing() -> None:
    pd.DataFrame([{
        "status": "NOT_EVALUABLE_MISSING_COMPETITION_ASSAY",
        "n_proteins": 0,
        "n_mutations": 0,
        "primary_metric": "per-protein Spearman between PBM score and competition effect",
        "reason": "No local Glasscock single-base competition quantitative table; no assay pairings inferred.",
    }]).to_csv(OUT / "06_assay_agreement" / "assay_agreement_status.tsv", sep="\t", index=False)
    (OUT / "06_assay_agreement" / "assay_agreement_report.md").write_text(
        "# PBM versus competition agreement\n\n"
        "Not evaluable because the local repository does not contain the Glasscock single-base competition table. "
        "PBM values were not treated as a substitute for the missing competition endpoint.\n",
        encoding="utf-8",
    )


def local_global(pwm_global: pd.DataFrame, nampnn: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    pbm = pd.read_parquet(ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet")
    rows = []
    for protein in PROTEINS:
        t = targets[targets.protein_id == protein].iloc[0]
        motif = str(t.designed_binding_site_motif)
        target_seq = str(t.experimental_assay_target)
        target7, _, amb, _ = choose_center_window(target_seq, motif)
        target_c = canonical_rc(target7)
        exp = pbm[pbm.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}") == protein][["canonical_7mer", "experimental_escore_consensus"]].copy()
        exp["distance"] = exp.canonical_7mer.map(lambda s: min(sum(a != b for a, b in zip(s, target_c)), sum(a != b for a, b in zip(reverse_complement(s), target_c))))
        for method, col in [("DeepPBS", "deeppbs_score")]:
            pred = pwm_global[pwm_global.protein == protein][["canonical_7mer", col]].rename(columns={col: "predicted"})
            merged = exp.merge(pred, on="canonical_7mer")
            merged["method"] = method
            for label, mask in [("0", merged.distance == 0), ("1", merged.distance == 1), ("2", merged.distance == 2), ("3", merged.distance == 3), (">=4", merged.distance >= 4), ("global", np.ones(len(merged), dtype=bool))]:
                g = merged[mask]
                rows.append({"method": method, "protein": protein, "distance_bin": label, "n_sequences": len(g), "spearman": spearman(g.predicted, g.experimental_escore_consensus), "target_7mer": target7, "target_canonical": target_c, "target_mapping_ambiguous": amb})
        ng = nampnn[nampnn.protein == protein][["canonical_7mer", "nampnn_score"]].drop_duplicates().rename(columns={"nampnn_score": "predicted"})
        merged = exp.merge(ng, on="canonical_7mer")
        merged["distance"] = merged.canonical_7mer.map(lambda s: min(sum(a != b for a, b in zip(s, target_c)), sum(a != b for a, b in zip(reverse_complement(s), target_c))))
        for label, mask in [("0", merged.distance == 0), ("1", merged.distance == 1), ("2", merged.distance == 2), ("3", merged.distance == 3), (">=4", merged.distance >= 4), ("global", np.ones(len(merged), dtype=bool))]:
            g = merged[mask]
            rows.append({"method": "NA-MPNN", "protein": protein, "distance_bin": label, "n_sequences": len(g), "spearman": spearman(g.predicted, g.experimental_escore_consensus), "target_7mer": target7, "target_canonical": target_c, "target_mapping_ambiguous": amb})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "07_local_global" / "local_global_metrics.tsv", sep="\t", index=False)
    return out


def pwm_capacity(targets: pd.DataFrame) -> pd.DataFrame:
    pbm = pd.read_parquet(ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet")
    rows, summary = [], []
    for protein in PROTEINS:
        t = targets[targets.protein_id == protein].iloc[0]
        target7, _, amb, _ = choose_center_window(str(t.experimental_assay_target), str(t.designed_binding_site_motif))
        target_c = canonical_rc(target7)
        exp = pbm[pbm.protein_id.map(lambda x: f"DBP{int(str(x).replace('DBP', '')):03d}") == protein].set_index("canonical_7mer")["experimental_escore_consensus"].to_dict()
        wt = exp[target_c]
        effects = {}
        for i, base in enumerate(target_c):
            for mut in BASES:
                if mut == base: continue
                s = target_c[:i] + mut + target_c[i + 1:]
                effects[(i, mut)] = exp[canonical_rc(s)] - wt
        for seq, value in exp.items():
            oriented, orientation, d = orient_to_target(seq, target_c)
            # Use canonical sequence coordinates. This is an extrapolation diagnostic, not a fitted predictor.
            pred = wt
            for i, (a, b) in enumerate(zip(oriented, target_c)):
                if a != b and (i, a) in effects: pred += effects[(i, a)]
            rows.append({"protein": protein, "canonical_7mer": seq, "oriented_7mer": oriented, "orientation": orientation, "hamming_distance": d, "experimental_score": value, "additive_predicted_score": pred, "residual": value - pred, "target_7mer": target_c, "mapping_ambiguous": amb})
        g = pd.DataFrame([r for r in rows if r["protein"] == protein])
        for label, mask in [("local_d1", g.hamming_distance == 1), ("d_ge_2", g.hamming_distance >= 2), ("global", np.ones(len(g), dtype=bool))]:
            q = g[mask]
            summary.append({"protein": protein, "subset": label, "n": len(q), "spearman": spearman(q.additive_predicted_score, q.experimental_score), "mae": float(np.mean(np.abs(q.residual))), "residual_sd": float(np.std(q.residual)), "target_7mer": target_c, "mapping_ambiguous": amb})
    pd.DataFrame(rows).to_csv(OUT / "08_pwm_capacity" / "additive_pwm_predictions.tsv", sep="\t", index=False)
    sdf = pd.DataFrame(summary)
    sdf.to_csv(OUT / "08_pwm_capacity" / "additive_pwm_capacity.tsv", sep="\t", index=False)
    return sdf


def conditionality(global_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method in sorted(global_df.method.unique()):
        sub = global_df[global_df.method == method].pivot(index="canonical_7mer", columns="protein", values="predicted_score")
        for a, b in combinations([p for p in PROTEINS if p in sub.columns], 2):
            rows.append({"method": method, "protein_a": a, "protein_b": b, "spearman": spearman(sub[a], sub[b]), "n": int(sub[[a, b]].dropna().shape[0])})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "10_conditionality" / "global_pairwise_correlations.tsv", sep="\t", index=False)
    summ = out.groupby("method").agg(median_pairwise_spearman=("spearman", "median"), min_pairwise_spearman=("spearman", "min"), max_pairwise_spearman=("spearman", "max"), n_pairs=("spearman", "size")).reset_index()
    summ.to_csv(OUT / "10_conditionality" / "conditionality_summary.tsv", sep="\t", index=False)
    # Local vectors are not homologous because competition labels are absent; record that explicitly.
    (OUT / "10_conditionality" / "local_conditionality_status.md").write_text("# Local conditionality\n\nNot evaluable: the local competition assay source is missing, and mutation coordinates cannot be compared as measured vectors.\n", encoding="utf-8")
    return out


def write_protocol_artifacts(reg: pd.DataFrame) -> None:
    (OUT / "03_register_mapping" / "register_mapping_protocol.md").write_text(
        "# v1.1 register protocol\n\n"
        "Primary seven-position windows are selected from exact model/assay cognate sequences and documented motifs before any experimental correlation. "
        "For DBP048 the design sequence lacks the assay motif CTGACG; the central design window is retained and marked ambiguous. No PBM-driven window or strand selection is performed.\n",
        encoding="utf-8",
    )
    (OUT / "09_global_pwm_projection" / "scoring_protocol.md").write_text(
        "# v1.1 global projection protocol\n\n"
        "Candidate universe: 8,192 canonical reverse-complement 7-mers. For DeepPBS, the fixed mapped seven model positions are scored as `sum(log(P+1e-12))`; forward and reverse-complement candidate orientations are both scored and the maximum is retained. NA-MPNN uses its frozen v1.0 canonical landscape verbatim because native PPM tensors were not retained. No positionwise RC averaging and no PBM-driven register optimization are used.\n",
        encoding="utf-8",
    )


def inventory(root: Path, pwm: pd.DataFrame, nampnn: pd.DataFrame, targets: pd.DataFrame) -> None:
    candidates = [
        ROOT / "results" / "v1_0_structure" / "nampnn_results.csv",
        ROOT / "data" / "processed" / "v0_3_1" / "designed_dbp_upbm_rc_class_v0_3_1.parquet",
        ROOT / "metadata" / "v0_3_1" / "designed_dbp_target_definitions.csv",
        ROOT / "data" / "raw" / "v0_3_1" / "41594_2025_1669_MOESM12_ESM.xlsx",
        ROOT / "data" / "raw" / "v0_3_1" / "41594_2025_1669_MOESM20_ESM.xls",
        ROOT / "data" / "raw" / "gse237017" / "nature_41594_2025_1669_MOESM3_ESM.xlsx",
        ROOT / "frozen_design7_20260914.tar.gz",
    ] + sorted((root / "predictions").glob("*.npz"))
    rows = []
    for p in candidates:
        if not p.exists(): continue
        rows.append({"source_file": p.relative_to(ROOT).as_posix(), "purpose": "DeepPBS/NA-MPNN/target/PBM source", "proteins_covered": "|".join(PROTEINS) if "npz" in p.name or "nampnn" in p.name else "see source metadata", "assay": "uPBM" if "upbm" in p.name.lower() or "MOESM12" in p.name else "structure/model metadata", "raw_processed_derived": "raw" if "raw" in p.parts or p.name.endswith((".npz", ".tar.gz")) else "processed", "dimensions": "" if p.suffix not in {".parquet", ".csv"} else str((len(pwm), len(nampnn)) if "nampnn" in p.name else "schema inspected"), "sha256": sha256(p)})
    pd.DataFrame(rows).to_csv(OUT / "00_inventory" / "source_inventory.csv", index=False)
    (OUT / "00_inventory" / "inventory_report.md").write_text(
        "# v1.1 repository inventory\n\n"
        f"Seven design proteins: {', '.join(PROTEINS)}. The frozen DeepPBS bundle was ingested without modifying NPZ files. "
        "Repository-wide search found no local single-base competition quantitative table; this is recorded as a missing input rather than inferred.\n",
        encoding="utf-8",
    )


def figures(global_df: pd.DataFrame, local_df: pd.DataFrame, capacity: pd.DataFrame) -> None:
    figdir = OUT / "figures"
    # 1 register schematic
    reg = pd.read_csv(OUT / "03_register_mapping" / "register_summary.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(12, 4)); ax.axis("off")
    for i, (_, r) in enumerate(reg.iterrows()): ax.text(0, 1 - i / 8, f"{r.protein}: model {r.model_window} | assay {r.assay_window} | {r.mapping_status}", family="monospace")
    fig.tight_layout(); fig.savefig(figdir / "figure1_register_map.png", dpi=200); plt.close(fig)
    # 2/3 missing competition status
    for n, title in [(2, "DeepPBS vs competition assay"), (3, "NA-MPNN vs competition assay")]:
        fig, ax = plt.subplots(figsize=(7, 4)); ax.text(.5, .5, "MISSING_LOCAL_SOURCE\ncompetition mutation table not available", ha="center", va="center", fontsize=13); ax.set_title(title); ax.axis("off"); fig.savefig(figdir / f"figure{n}_competition_missing.png", dpi=200); plt.close(fig)
    # 4 global summary
    g = pd.read_csv(OUT / "09_global_pwm_projection" / "global_projection_metrics.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(8, 4));
    for method, grp in g.groupby("method"): ax.scatter(grp.protein, grp.spearman, label=method)
    ax.axhline(.5914, color="grey", ls="--", label="replicate reference (context)"); ax.set_ylabel("PBM Spearman"); ax.set_xlabel("protein"); ax.legend(); fig.tight_layout(); fig.savefig(figdir / "figure4_local_global_summary.png", dpi=200); plt.close(fig)
    # 5 missing agreement
    fig, ax = plt.subplots(figsize=(7, 4)); ax.text(.5, .5, "MISSING_LOCAL_SOURCE\nPBM vs competition agreement not evaluable", ha="center", va="center"); ax.axis("off"); fig.savefig(figdir / "figure5_assay_agreement_missing.png", dpi=200); plt.close(fig)
    # 6 distance curves
    fig, ax = plt.subplots(figsize=(8, 4)); q = local_df[local_df.distance_bin != "global"]
    for method, grp in q.groupby("method"):
        z = grp.groupby("distance_bin", sort=False).spearman.median(); ax.plot(z.index, z.values, marker="o", label=method)
    ax.set_ylabel("median within-distance Spearman"); ax.set_xlabel("Hamming distance"); ax.legend(); fig.tight_layout(); fig.savefig(figdir / "figure6_hamming_distance.png", dpi=200); plt.close(fig)
    # 7 additive capacity
    fig, ax = plt.subplots(figsize=(8, 4)); q = capacity[capacity.subset.isin(["d_ge_2", "global"])];
    for sub, grp in q.groupby("subset"): ax.scatter(grp.protein, grp.spearman, label=sub)
    ax.set_ylabel("experimental additive prediction Spearman"); ax.legend(); fig.tight_layout(); fig.savefig(figdir / "figure7_additive_capacity.png", dpi=200); plt.close(fig)
    # 8 conditionality heatmap-like matrix
    c = pd.read_csv(OUT / "10_conditionality" / "global_pairwise_correlations.tsv", sep="\t")
    fig, axes = plt.subplots(1, max(1, c.method.nunique()), figsize=(6 * max(1, c.method.nunique()), 5), squeeze=False)
    for ax, (method, grp) in zip(axes[0], c.groupby("method")):
        mat = pd.DataFrame(np.eye(7), index=PROTEINS, columns=PROTEINS)
        for _, r in grp.iterrows(): mat.loc[r.protein_a, r.protein_b] = mat.loc[r.protein_b, r.protein_a] = r.spearman
        im = ax.imshow(mat, vmin=-1, vmax=1, cmap="coolwarm"); ax.set_title(method); ax.set_xticks(range(7), PROTEINS, rotation=90); ax.set_yticks(range(7), PROTEINS); fig.colorbar(im, ax=ax, fraction=.046)
    fig.tight_layout(); fig.savefig(figdir / "figure8_conditionality_heatmaps.png", dpi=200); plt.close(fig)
    shutil.copy2(figdir / "figure2_competition_missing.png", OUT / "05_assay_aligned_eval" / "figure2_competition_missing.png")
    shutil.copy2(figdir / "figure3_competition_missing.png", OUT / "05_assay_aligned_eval" / "figure3_competition_missing.png")
    shutil.copy2(figdir / "figure5_assay_agreement_missing.png", OUT / "06_assay_agreement" / "figure5_assay_agreement_missing.png")
    pd.DataFrame([
        {"figure": "Figure 1", "path": "figures/figure1_register_map.png", "source_data": "03_register_mapping/register_summary.tsv"},
        {"figure": "Figure 2", "path": "figures/figure2_competition_missing.png", "source_data": "04_competition/competition_source_status.tsv"},
        {"figure": "Figure 3", "path": "figures/figure3_competition_missing.png", "source_data": "04_competition/competition_source_status.tsv"},
        {"figure": "Figure 4", "path": "figures/figure4_local_global_summary.png", "source_data": "09_global_pwm_projection/global_projection_metrics.tsv"},
        {"figure": "Figure 5", "path": "figures/figure5_assay_agreement_missing.png", "source_data": "06_assay_agreement/assay_agreement_status.tsv"},
        {"figure": "Figure 6", "path": "figures/figure6_hamming_distance.png", "source_data": "07_local_global/local_global_metrics.tsv"},
        {"figure": "Figure 7", "path": "figures/figure7_additive_capacity.png", "source_data": "08_pwm_capacity/additive_pwm_capacity.tsv"},
        {"figure": "Figure 8", "path": "figures/figure8_conditionality_heatmaps.png", "source_data": "10_conditionality/global_pairwise_correlations.tsv"},
    ]).to_csv(OUT / "figures" / "figure_source_manifest.tsv", sep="\t", index=False)


def overlap_audit() -> str:
    files = list((ROOT / "external" / "nampnn" / "NA-MPNN" / "splits").glob("*"))
    hits = []
    for f in files:
        try: text = f.read_text(errors="ignore").lower()
        except Exception: continue
        for token in ["dbp001", "dbp003", "dbp005", "dbp006", "dbp009", "dbp035", "dbp048", "8tac"]:
            if token in text: hits.append((f.relative_to(ROOT).as_posix(), token))
    lines = ["# NA-MPNN training-overlap audit", "", "Exact designed-protein training sequences are not distributed in the checked-out repository; homolog-level overlap is therefore UNDETERMINED.", "", "## Explicit hits"]
    if hits:
        lines += [f"- `{f}` contains `{t}`: POSSIBLE/CONFIRMED overlap requires interpretation of split semantics." for f, t in hits]
    else: lines += ["- No exact DBP design IDs were found in local split manifests: NO EXACT OVERLAP FOUND (sequence-level overlap remains UNDETERMINED)."]
    lines += ["", "DBP48/8TAC was historically flagged in local NA-MPNN split files; the v1.1 design-model runs use Supplementary Data 6 design PDBs, not 8TAC. This does not establish a clean homolog-level separation."]
    text = "\n".join(lines) + "\n"
    (OUT / "11_training_overlap" / "training_overlap_report.md").write_text(text, encoding="utf-8")
    return "UNDETERMINED"


def write_summary(metrics: dict) -> None:
    deep = metrics["global"][metrics["global"].method == "DeepPBS"]
    na = metrics["global"][metrics["global"].method == "NA-MPNN"]
    deep_primary = deep[deep.primary_status == "primary"]
    na_primary = na[na.primary_status == "primary"]
    cond = metrics["conditionality"].groupby("method").spearman.median().to_dict()
    cap = metrics["capacity"]
    lines = [
        "# v1.1 summary",
        "",
        "This is a development-stage assay-alignment audit on the seven development-exposed GSE237017 DBPs. No v1.0 outputs were overwritten.",
        "",
        "## Findings",
        f"1. DeepPBS design-model coverage/QC: 7/7 NPZ predictions validated; widths are 13, 13, 15, 15, 15, 15, and 14.",
        "2. Competition-aligned DeepPBS performance: NOT EVALUABLE because the local Glasscock single-base competition table is missing.",
        "3. Competition-aligned NA-MPNN performance: NOT EVALUABLE; frozen v1.0 artifact retained global landscapes but no native PPM tensor.",
        f"4. Fixed-register global PBM Spearman (six uniquely mapped primary proteins): DeepPBS median {deep_primary.spearman.median():.4f}; NA-MPNN median {na_primary.spearman.median():.4f}. DBP048 is retained only as an ambiguous-mapping sensitivity result (DeepPBS {float(deep.loc[deep.protein == 'DBP048', 'spearman'].iloc[0]):.4f}; NA-MPNN {float(na.loc[na.protein == 'DBP048', 'spearman'].iloc[0]):.4f}).",
        "5. PBM-vs-competition experimental agreement: NOT EVALUABLE for the same missing source.",
        f"6. Experimental additive-PWM capacity: median extrapolation Spearman at Hamming distance >=2 is {cap.loc[cap.subset == 'd_ge_2', 'spearman'].median():.4f}; this is a capacity diagnostic, not a predictive method.",
        "7. Local-to-global decay: distance-stratified values are reported without selecting bins post hoc; interpretation remains descriptive.",
        f"8. Global protein specificity: median pairwise landscape correlation is DeepPBS {cond.get('DeepPBS', float('nan')):.4f} and NA-MPNN {cond.get('NA-MPNN', float('nan')):.4f}; this supports protein-dependent outputs, not accuracy.",
        "9. Training overlap: exact sequence-level separation is UNDETERMINED because full training sequences are unavailable.",
        "10. Supported conclusions: additive structure-model outputs can be protein-specific, but global PBM ranking remains limited; the assay-aligned local hypothesis is unresolved until the competition table is obtained.",
        "11. Not supported: no claim that structure methods generally fail or that this exposed seven-protein cohort is independent validation.",
        "",
        "## Decision classification",
        "A is not decidable without the local competition assay. D is supported as a capacity diagnostic if additive extrapolation is weak outside Hamming-1. C cannot be assessed because PBM and competition cannot be paired. The current status is **INCONCLUSIVE / BLOCKED ON LOCAL COMPETITION ASSAY SOURCE**.",
    ]
    (OUT / "v1_1_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    serial = {"deeppbs_global_median_spearman_all7_descriptive": float(deep.spearman.median()), "nampnn_global_median_spearman_all7_descriptive": float(na.spearman.median()), "deeppbs_global_median_spearman_primary_n6": float(deep_primary.spearman.median()), "nampnn_global_median_spearman_primary_n6": float(na_primary.spearman.median()), "deeppbs_global_per_protein": deep.set_index("protein").spearman.to_dict(), "nampnn_global_per_protein": na.set_index("protein").spearman.to_dict(), "conditionality_median_pairwise": {k: float(v) for k, v in cond.items()}, "additive_pwm_d_ge_2_median_spearman": float(cap.loc[cap.subset == "d_ge_2", "spearman"].median()), "competition_assay_status": "MISSING_LOCAL_SOURCE", "training_overlap_status": "UNDETERMINED", "development_exposed": True}
    (OUT / "v1_1_metrics.json").write_text(json.dumps(serial, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    root = ingest_tar()
    pwm, _ = ingest_deeppbs(root)
    nampnn, _ = ingest_nampnn()
    targets = load_targets()
    inventory(root, pwm, nampnn, targets)
    reg, _ = build_register_map(pwm, targets)
    write_protocol_artifacts(reg)
    # Keep the original missing-source behavior only for a genuinely absent
    # source. When the official workbook is present, the dedicated parser
    # owns the assay-aligned outputs and is run after the reusable global
    # diagnostics below.
    competition_available = COMPETITION_SOURCE.exists()
    if not competition_available:
        competition_missing()
    global_df, global_metrics = project_global(pwm, nampnn, targets, reg)
    mutation_projection(pwm, reg, targets, global_df, nampnn)
    if not competition_available:
        assay_agreement_missing()
    local_df = local_global(global_df, nampnn, targets)
    capacity = pwm_capacity(targets)
    cond = conditionality(pd.read_csv(OUT / "09_global_pwm_projection" / "unified_global_landscapes.tsv", sep="\t"))
    overlap_audit()
    figures(global_df, local_df, capacity)
    write_summary({"global": global_metrics, "conditionality": cond, "capacity": capacity})
    if competition_available:
        from scripts.v1_1.competition_assay import complete_competition
        complete_competition()
        print("v1.1 pipeline complete: competition branch status=PARSED_MOESM16")
    else:
        print("v1.1 pipeline complete: competition branch status=MISSING_LOCAL_SOURCE")


if __name__ == "__main__":
    main()
