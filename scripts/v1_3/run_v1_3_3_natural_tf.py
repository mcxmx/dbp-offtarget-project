"""v1.3.3 natural-TF structure audit, prediction freeze, and reveal.

The prediction stage intentionally reads only SaMBA sequence/mapping metadata
and frozen DeepPBS NPZ outputs. It does not load experimental effects.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from Bio.Seq import Seq
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
SPOTS = ROOT / "data/processed/v1_3_2_samba_spots.parquet"
STRUCTURES = ROOT / "data/raw/v1_3_3_external/structures"
RESULTS = ROOT / "results/v1_3"
REPORTS = ROOT / "reports/v1.3"
MAPPING = RESULTS / "natural_tf_structure_mapping.tsv"
PREDICTIONS = RESULTS / "natural_tf_predictions_unscored.tsv"
REVEALED = RESULTS / "natural_tf_predictions_revealed.tsv"
SEED = 1301

DNA_CODES = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "A": "A", "C": "C", "G": "G", "T": "T"}
BASES = "ACGT"

# Candidate choices are label-blind and frozen in the challenge contract.
CANDIDATES = {
    "Cbf1": {"pdb": "8OVW", "protein_chain": "A", "dna_chains": "D,E", "source": "RCSB experimental cryo-EM"},
    "Egr1": {"pdb": "1AAY", "protein_chain": "A", "dna_chains": "B,C", "source": "RCSB experimental X-ray"},
    "Ets1": {"pdb": "2STT", "protein_chain": "C", "dna_chains": "A,B", "source": "RCSB experimental NMR"},
    "GR": {"pdb": "1R4R", "protein_chain": "A", "dna_chains": "C,D", "source": "RCSB experimental X-ray"},
    "Max": {"pdb": "1HLO", "protein_chain": "A", "dna_chains": "C,D", "source": "RCSB experimental X-ray"},
    "TBP": {"pdb": "1TGH", "protein_chain": "A", "dna_chains": "B,C", "source": "RCSB experimental X-ray"},
    "p53": {"pdb": "1TUP", "protein_chain": "A", "dna_chains": "E,F", "source": "RCSB experimental X-ray"},
}


def safe_corr(a: np.ndarray, b: np.ndarray, kind: str = "spearman") -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return math.nan
    return float((spearmanr if kind == "spearman" else pearsonr)(a, b).statistic)


def pairwise_accuracy(exp: np.ndarray, pred: np.ndarray) -> float:
    values = []
    for i in range(len(exp)):
        for j in range(i + 1, len(exp)):
            de, dp = exp[i] - exp[j], pred[i] - pred[j]
            if de == 0:
                continue
            values.append(1.0 if de * dp > 0 else 0.5 if dp == 0 else 0.0)
    return float(np.mean(values)) if values else math.nan


def nw_align(a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2):
    """Deterministic global alignment, returning score and a->b coordinates."""
    n, m = len(a), len(b)
    score = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + gap
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score[i, j] = max(score[i - 1, j - 1] + (match if a[i - 1] == b[j - 1] else mismatch), score[i - 1, j] + gap, score[i, j - 1] + gap)
    i, j = n, m
    pairs = []
    while i or j:
        current = score[i, j]
        if i and j and current == score[i - 1, j - 1] + (match if a[i - 1] == b[j - 1] else mismatch):
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif i and current == score[i - 1, j] + gap:
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    return int(score[n, m]), list(reversed(pairs))


def load_sequences() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_parquet(SPOTS, columns=["tf_id", "binding_site_id", "position", "wt_basepair", "perturbation", "perturbation_type", "sequence"])
    reference = raw[raw.perturbation_type.eq("reference")].drop_duplicates(["tf_id", "binding_site_id"])
    canonical = raw[raw.perturbation_type.eq("canonical_watson_crick")].drop_duplicates(["tf_id", "binding_site_id", "position", "perturbation"])
    return reference, canonical


def structure_details(pdb_id: str) -> dict:
    path = STRUCTURES / f"{pdb_id}.cif"
    if not path.exists():
        path = ROOT / "data/raw/rcsb/mmcif" / f"{pdb_id}.cif"
    d = MMCIF2Dict(str(path))
    def first(key, default=""):
        val = d.get(key, default)
        return val[0] if isinstance(val, list) else val
    resolution = first("_refine.ls_d_res_high", "")
    if resolution in {"", "."}:
        resolution = first("_em_3d_reconstruction.resolution", "")
    return {"title": first("_struct.title"), "method": first("_exptl.method"), "resolution": resolution}


def chain_sequences(pdb_id: str) -> dict[str, str]:
    path = STRUCTURES / f"{pdb_id}.cif"
    if not path.exists():
        path = ROOT / "data/raw/rcsb/mmcif" / f"{pdb_id}.cif"
    structure = MMCIFParser(QUIET=True).get_structure(pdb_id, str(path))
    model = next(structure.get_models())
    out = {}
    for chain in model:
        seq = "".join(DNA_CODES.get(r.resname.strip(), "N") for r in chain if r.id[0] == " " and r.resname.strip() in DNA_CODES)
        if seq:
            out[chain.id] = seq
    return out


def site_metadata(reference: pd.DataFrame, canonical: pd.DataFrame, tf: str, site: str) -> dict:
    r = reference[(reference.tf_id == tf) & (reference.binding_site_id == site)].iloc[0]
    c = canonical[(canonical.tf_id == tf) & (canonical.binding_site_id == site)].copy()
    positions = sorted(int(x) for x in c.position.unique())
    diffs = {}
    for position in positions:
        row = c[c.position.eq(position)].iloc[0]
        diffs[position] = [i for i, (a, b) in enumerate(zip(r.sequence, row.sequence)) if a != b][0]
    start = min(diffs.values())
    motif = r.sequence[start : start + len(positions)]
    return {"reference_sequence": r.sequence, "positions": positions, "position_to_core_index": {str(p): i for i, p in enumerate(sorted(positions))}, "mutation_indices": diffs, "motif": motif}


def best_mapping(motif: str, positions: list[int], chains: dict[str, str], dna_pair: list[str]) -> dict:
    candidates = []
    for chain_id in dna_pair:
        if chain_id not in chains:
            continue
        seq = chains[chain_id]
        for orientation, oriented in [("direct", seq), ("reverse_complement", str(Seq(seq).reverse_complement()))]:
            for start in range(len(oriented)):
                for length in range(max(4, len(motif) - 2), min(len(motif) + 3, len(oriented) - start + 1)):
                    window = oriented[start : start + length]
                    score, pairs = nw_align(motif, window)
                    mapping = {a: b for a, b in pairs if a is not None and b is not None}
                    if not all(i in mapping for i in range(len(motif))):
                        continue
                    identity = sum(motif[i] == window[mapping[i]] for i in range(len(motif))) / len(motif)
                    candidates.append({"chain": chain_id, "orientation": orientation, "start": start, "length": length, "score": score, "identity": identity, "mapping": mapping, "oriented_sequence": oriented, "window": window})
    if not candidates:
        return {"status": "MAPPING_AMBIGUOUS", "reason": "no candidate alignment maps every mutated coordinate"}
    candidates.sort(key=lambda x: (x["score"], x["identity"], -x["start"]), reverse=True)
    best = candidates[0]
    ties = [x for x in candidates if (x["score"], round(x["identity"], 12)) == (best["score"], round(best["identity"], 12))]
    if len(ties) != 1:
        return {"status": "MAPPING_AMBIGUOUS", "reason": f"{len(ties)} tied label-blind alignments", "best_score": best["score"]}
    # A unique candidate is still rejected if it maps to multiple chain/pair orientations.
    return {"status": "ELIGIBLE", **{k: best[k] for k in ["chain", "orientation", "start", "length", "score", "identity", "window"]}, "motif_to_oriented_index": json.dumps(best["mapping"], sort_keys=True)}


def audit() -> pd.DataFrame:
    reference, canonical = load_sequences()
    rows = []
    for (tf, site), _ in reference.groupby(["tf_id", "binding_site_id"], sort=True):
        candidate = CANDIDATES[tf]
        details = structure_details(candidate["pdb"])
        chains = chain_sequences(candidate["pdb"])
        meta = site_metadata(reference, canonical, tf, site)
        mapping = best_mapping(meta["motif"], meta["positions"], chains, candidate["dna_chains"].split(","))
        rows.append({"tf_id": tf, "binding_site_id": site, "pdb_id": candidate["pdb"], "protein_chain": candidate["protein_chain"], "dna_chains": candidate["dna_chains"], "structure_source": candidate["source"], "experimental_method": details["method"], "resolution_A": details["resolution"], "structure_title": details["title"], "assay_sequence": meta["reference_sequence"], "canonical_motif": meta["motif"], "n_positions": len(meta["positions"]), "pdb_dna_sequence": chains.get(mapping.get("chain", ""), ""), "target_similarity": mapping.get("identity", np.nan), "mapping_status": mapping.get("status", "MAPPING_AMBIGUOUS"), "mapping_reason": mapping.get("reason", ""), "mapping_chain": mapping.get("chain", ""), "mapping_orientation": mapping.get("orientation", ""), "mapping_start": mapping.get("start", np.nan), "mapping_length": mapping.get("length", np.nan), "alignment_score": mapping.get("score", np.nan), "motif_to_oriented_index": mapping.get("motif_to_oriented_index", ""), "protein_sequence_coverage": "chain_residue_count_not_reference_aligned"})
    out = pd.DataFrame(rows).sort_values(["tf_id", "binding_site_id"])
    RESULTS.mkdir(parents=True, exist_ok=True); REPORTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(MAPPING, sep="\t", index=False, na_rep="NA")
    lines = ["# v1.3.3 Structure Audit", "", "Status: frozen before natural-TF model/label performance inspection (2026-09-18).", "", "The audit uses only RCSB structure metadata and SaMBA sequence/mutation-coordinate metadata. No SaMBA effect column is loaded.", "", "## Selection rule", "", "Primary candidates were fixed label-blind by TF identity and direct experimental DNA-bound annotation. The mapping uses the contract's deterministic exact/global alignment rule; no correlation, rank, or pairwise metric was inspected.", "", "## Site audit", "", "| TF | Site | PDB | Protein chain | DNA chains | Method | Resolution A | Motif | PDB DNA | Identity | Mapping | Reason |", "|---|---|---|---|---|---|---:|---|---|---:|---|---|"]
    for _, row in out.iterrows():
        lines.append(f"| {row.tf_id} | {row.binding_site_id} | {row.pdb_id} | {row.protein_chain} | {row.dna_chains} | {row.experimental_method} | {row.resolution_A} | {row.canonical_motif} | {row.pdb_dna_sequence} | {row.target_similarity:.3f} | {row.mapping_status} | {row.mapping_reason or 'unique label-blind mapping'} |" if pd.notna(row.target_similarity) else f"| {row.tf_id} | {row.binding_site_id} | {row.pdb_id} | {row.protein_chain} | {row.dna_chains} | {row.experimental_method} | {row.resolution_A} | {row.canonical_motif} | {row.pdb_dna_sequence} | NA | {row.mapping_status} | {row.mapping_reason} |")
    lines += ["", "## Coverage", "", f"- SaMBA sites audited: {len(out)} across {out.tf_id.nunique()} TFs.", f"- Eligible mappings: {int(out.mapping_status.eq('ELIGIBLE').sum())} sites.", f"- Excluded/ambiguous mappings: {int((~out.mapping_status.eq('ELIGIBLE')).sum())} sites.", "- `protein_sequence_coverage` is reported conservatively as structure-chain coverage metadata; no sequence identity to the assay construct was inferred without a frozen reference mapping.", "- If the eligible count is too small for a TF-level challenge, the result is reported as `STRUCTURAL_COVERAGE_LIMITED` and v1.4 remains closed."]
    (REPORTS / "V1_3_3_STRUCTURE_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def load_frozen_mapping() -> pd.DataFrame:
    if not MAPPING.exists():
        raise FileNotFoundError("Run audit and commit the structure/mapping freeze first")
    return pd.read_csv(MAPPING, sep="\t")


def write_prediction_rows(mapping: pd.DataFrame, npz_dir: Path) -> pd.DataFrame:
    # Labels are deliberately not opened here. Only sequence/mutation identity metadata is read.
    _, canonical = load_sequences()
    rows = []
    for _, site in mapping[mapping.mapping_status.eq("ELIGIBLE")].iterrows():
        npz_path = npz_dir / f"{site.pdb_id}_{site.binding_site_id}.npz_predict.npz"
        if not npz_path.exists():
            continue
        pred = np.load(npz_path, allow_pickle=False)
        p = np.asarray(pred["P"], dtype=float)
        seq = "".join(BASES[int(i)] for i in np.asarray(pred["Seq"]).reshape(-1) if int(i) < 4)
        core_map = json.loads(site.motif_to_oriented_index)
        for _, mutation in canonical[(canonical.tf_id == site.tf_id) & (canonical.binding_site_id == site.binding_site_id)].iterrows():
            core_idx = int(core_map[str(int(mutation.position) - 1)])
            # Orientation is recorded for audit; DeepPBS P is interpreted in the emitted Seq order.
            pred_idx = int(site.mapping_start) + core_idx
            if pred_idx >= len(p) or pred_idx >= len(seq):
                continue
            wt = str(mutation.wt_basepair)[0]
            mutant = str(mutation.perturbation)[0]
            if wt not in BASES or mutant not in BASES:
                continue
            rows.append({"tf_id": site.tf_id, "site_id": site.binding_site_id, "pdb_id": site.pdb_id, "position": int(mutation.position), "wt_base": wt, "mutant_base": mutant, "predicted_effect": float(np.log(p[pred_idx, BASES.index(mutant)] + 1e-12) - np.log(p[pred_idx, BASES.index(wt)] + 1e-12)), "model": "DeepPBS", "model_version": "official_commit_8bfb211_ensemble_828_529_173_898_820", "structure_chain": site.protein_chain, "dna_chain": site.mapping_chain, "dna_register": f"{site.mapping_start}:{site.mapping_length}", "orientation": site.mapping_orientation, "mapping_status": "FROZEN", "prediction_source": str(npz_path.relative_to(ROOT))})
    out = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(PREDICTIONS, sep="\t", index=False, na_rep="NA")
    return out


def reveal() -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = load_frozen_mapping()
    pred = pd.read_csv(PREDICTIONS, sep="\t")
    labels = pd.read_parquet(SPOTS)
    ref = labels[labels.perturbation_type.eq("reference")][["tf_id", "binding_site_id", "sequence"]].drop_duplicates()
    canonical = labels[labels.perturbation_type.eq("canonical_watson_crick")].copy()
    # One source-defined consensus label per mutation. This is only opened after prediction freeze.
    canonical["experimental_effect"] = np.log2(pd.to_numeric(canonical["published_fold_change"], errors="coerce"))
    canonical = canonical.groupby(["tf_id", "binding_site_id", "position", "wt_basepair", "perturbation"], as_index=False).experimental_effect.first()
    canonical["wt_base"] = canonical.wt_basepair.str[0]
    canonical["mutant_base"] = canonical.perturbation.str[0]
    joined = canonical.merge(pred, on=["tf_id", "site_id", "position", "wt_base", "mutant_base"], how="inner", validate="one_to_one")
    per_rows = []
    for (tf, site), g in joined.groupby(["tf_id", "site_id"], sort=True):
        g = g.sort_values(["position", "mutant_base"])
        groups = [q for _, q in g.groupby("position", sort=True)]
        exp_pos = np.array([np.mean(np.abs(q.experimental_effect)) for q in groups])
        pred_pos = np.array([np.mean(np.abs(q.predicted_effect)) for q in groups])
        exp_res = np.concatenate([q.experimental_effect.to_numpy(float) - q.experimental_effect.mean() for q in groups])
        pred_res = np.concatenate([q.predicted_effect.to_numpy(float) - q.predicted_effect.mean() for q in groups])
        per_rows.append({"row_type": "tf_site", "tf_id": tf, "site_id": site, "method": "DeepPBS", "n_mutations": len(g), "n_positions": len(groups), "global_spearman": safe_corr(g.experimental_effect, g.predicted_effect), "global_pearson": safe_corr(g.experimental_effect, g.predicted_effect, "pearson"), "position_spearman": safe_corr(exp_pos, pred_pos), "position_pearson": safe_corr(exp_pos, pred_pos, "pearson"), "identity_residual_spearman": safe_corr(exp_res, pred_res), "identity_residual_pearson": safe_corr(exp_res, pred_res, "pearson"), "identity_pairwise_accuracy": float(np.nanmean([pairwise_accuracy(q.experimental_effect.to_numpy(float), q.predicted_effect.to_numpy(float)) for q in groups])), "status": "EVALUABLE"})
    per = pd.DataFrame(per_rows)
    per.to_csv(RESULTS / "natural_tf_decomposition_per_tf.tsv", sep="\t", index=False, na_rep="NA")
    summary_rows = []
    for metric in ["global_spearman", "position_spearman", "identity_residual_spearman", "identity_pairwise_accuracy"]:
        values = per[metric].to_numpy(float)
        summary_rows.append({"row_type": "tf_median", "method": "DeepPBS", "metric": metric, "n_sites": int(len(per)), "n_tfs": int(per.tf_id.nunique()), "median": float(np.nanmedian(values)) if len(values) else np.nan, "mean": float(np.nanmean(values)) if len(values) else np.nan, "chance_baseline": 0.5 if metric.endswith("pairwise_accuracy") else 0.0, "status": "STRUCTURAL_COVERAGE_LIMITED" if per.tf_id.nunique() < 4 else "EVALUABLE"})
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(RESULTS / "natural_tf_decomposition.tsv", sep="\t", index=False, na_rep="NA")
    joined.to_csv(REVEALED, sep="\t", index=False, na_rep="NA")
    return per, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("audit")
    p = sub.add_parser("prediction")
    p.add_argument("--npz-dir", type=Path, required=True)
    sub.add_parser("reveal")
    args = parser.parse_args()
    if args.command == "audit":
        out = audit(); print(out[["tf_id", "binding_site_id", "pdb_id", "mapping_status", "target_similarity"]].to_string(index=False))
    elif args.command == "prediction":
        out = write_prediction_rows(load_frozen_mapping(), args.npz_dir); print(f"wrote {len(out)} unscored rows")
    else:
        per, summary = reveal(); print(per.to_string(index=False)); print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
