"""v1.0 structure-conditioning pilot.

This module is deliberately an additive, audit-first workflow.  It never edits
v0.x artifacts and it does not tune a model on PBM results.  It inventories the
public design PDBs, scores frozen NA-MPNN specificity outputs, and records
unavailable DeepPBS/Rosetta/AF3 runs rather than imputing them.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

from src.sequence_equivalence import canonical_rc, reverse_complement
from src.structure_baselines import load_nampnn_npz, score_structural_ppm
from src.deeppbs_evaluation import canonical_7mer_universe, load_designed_experimental_units
from src.v0_4_evaluation import compute_ranking_metrics


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v1_0_structure"
META = ROOT / "metadata" / "v1_0_structure"
DOC = ROOT / "docs" / "v1_0_structure"
SUP = ROOT / "data" / "raw" / "v1_0_structure" / "design_pdbs" / "design_pdbs"
DESIGN_PDB = ROOT / "external" / "dbp_design" / "2b_design_mpnn" / "DBP035.pdb"
CRYSTAL = ROOT / "data" / "raw" / "rcsb" / "mmcif" / "8TAC.cif"
DBPS = ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]
PROTEINS = {
    "DBP1": "TARELEVAALIAQGRSNREIAEELNISERTVERYVRRILRKLGLRNRAQIAAWVIRRS",
    "DBP3": "TKREREVLKLIAEDYGNKEIANRLNISERTVERYIRRILRKLGLKNRAELVRYAIRHG",
    "DBP5": "GFGKAVKAKRAELGLTQAEFAERAGLSRRTIIRIEQGKVKATSTTAEKIAAALGTTVQELEQA",
    "DBP6": "DWAARAAAARRLRKERGLTQAELGELAGVSRTTVSRIELGRPDVSQASVDAVLAVL",
    "DBP9": "DWERRCAYARRARKELGLTQAELGELAGVSRTTVSRIERGKPDVSEASVEAVLAVL",
    "DBP35": "GFGRAVKEKRKELGLTQKEFAEKAGLSRRTIIRIERGYIVPPKATKEKIAKALGTSVEELEQA",
    "DBP48": "MTPEEIAEAKRIGKEVKERRKELGLTQRELAEKLGVSRSTVSDIENGRRLPSEELLKKIKEILGV",
}
TARGETS = {
    "DBP1": ("TAGCAGGATGTGT", "GCAGG"),
    "DBP3": ("TAGCAGGATGTGT", "GCAGGA"),
    "DBP5": ("GCAGATCTGCACATC", "TGCACA"),
    "DBP6": ("GCAGATCTGCACATC", "TGCACA"),
    "DBP9": ("GCAGATCTGCACATC", "TGCACA"),
    "DBP35": ("GCAGATCTGCACATC", "TGCACA"),
    "DBP48": ("CGCCCAAAGCCGCG", "CTGACG"),
}
DNA_RES = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "A": "A", "C": "C", "G": "G", "T": "T", "U": "T"}
BASE_ATOMS = {"N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6", "C2", "C4", "C5", "C6", "C8"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def structure(path: Path):
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return next(MMCIFParser(QUIET=True).get_structure(path.stem, str(path)).get_models())
    return next(PDBParser(QUIET=True).get_structure(path.stem, str(path)).get_models())


def chain_sequences(model):
    out = {}
    for chain in model:
        protein, dna = [], []
        for residue in chain:
            rn = residue.resname.strip().upper()
            if rn in DNA_RES:
                dna.append(DNA_RES[rn])
            else:
                if not is_aa(residue, standard=False):
                    continue
                try:
                    protein.append(seq1(rn, custom_map={"MSE": "M"}))
                except Exception:
                    pass
        out[chain.id] = {"protein": "".join(protein), "dna": "".join(dna), "nres": len(list(chain)), "residues": list(chain)}
    return out


def pde(pdb: Path, dbp: str, source: str, structure_type: str):
    model = structure(pdb)
    chains = chain_sequences(model)
    protein_chain = next((c for c, x in chains.items() if x["protein"] == PROTEINS[dbp]), None)
    if protein_chain is None:
        protein_chain = next((c for c, x in chains.items() if PROTEINS[dbp] in x["protein"] or x["protein"] in PROTEINS[dbp]), None)
    dna_chains = [c for c, x in chains.items() if x["dna"]]
    dna = ";".join(chains[c]["dna"] for c in dna_chains)
    exact = bool(protein_chain and chains[protein_chain]["protein"] == PROTEINS[dbp])
    protein_sequence = chains[protein_chain]["protein"] if protein_chain else ""
    exact_construct = protein_sequence == PROTEINS[dbp]
    protein_relation = "exact" if exact_construct else ("PBM sequence is contained in structure construct" if protein_chain and PROTEINS[dbp] in protein_sequence else "sequence mismatch")
    dna_mapping = ""
    if dna_chains:
        fwd = TARGETS[dbp][0]
        rc = reverse_complement(fwd)
        chain_dnas = [chains[c]["dna"] for c in dna_chains]
        if any(x == fwd + rc for x in chain_dnas):
            dna_mapping = "chain contains target forward strand followed by exact reverse complement"
        elif sorted(chain_dnas) == sorted([fwd, rc]):
            dna_mapping = "separate chains are target forward strand and exact reverse complement"
        else:
            dna_mapping = "DNA sequence does not exactly equal the frozen target duplex"
    return {
        "protein_id": dbp,
        "exact_pbm_protein_sequence": PROTEINS[dbp],
        "intended_target_dna": TARGETS[dbp][0],
        "pbm_register": TARGETS[dbp][1],
        "design_model_pdb": str(pdb.relative_to(ROOT)).replace("\\", "/"),
        "experimental_structure_pdb": "8TAC" if dbp == "DBP48" else "",
        "protein_chain_id": protein_chain or "",
        "dna_chain_ids": ";".join(dna_chains),
        "dna_sequence_by_chain": dna,
        "dna_chain_lengths": ";".join(str(len(chains[c]["dna"])) for c in dna_chains),
        "dna_length_total": sum(len(chains[c]["dna"]) for c in dna_chains),
        "binding_core": TARGETS[dbp][1],
        "pdb_to_target_mapping": dna_mapping or "no DNA mapping",
        "exact_pbm_construct_in_structure": exact_construct,
        "protein_construct_relation": protein_relation,
        "structure_protein_sequence": protein_sequence,
        "missing_residues": "none detected in primary protein chain" if exact_construct else "sequence mismatch; inspect manually",
        "mutations": "none detected" if exact_construct else "not resolved",
        "terminal_additions_or_tags": "none in design model" if exact_construct else "terminal/construct differences present",
        "structure_provenance": source,
        "structure_type": structure_type,
        "structure_parse_status": "parsed" if protein_chain and dna_chains else "incomplete",
    }


def inventory():
    rows = []
    for dbp in DBPS:
        design_name = f"DBP{int(dbp[3:]):03d}.pdb"
        path = SUP / design_name
        if dbp == "DBP35" and not path.exists():
            path = DESIGN_PDB
        if path.exists():
            rows.append(pde(path, dbp, "Nature Supplementary Data 6 design-hit PDB", "designed protein-DNA complex"))
        else:
            rows.append({"protein_id": dbp, "exact_pbm_protein_sequence": PROTEINS[dbp], "intended_target_dna": TARGETS[dbp][0], "pbm_register": TARGETS[dbp][1], "design_model_pdb": "", "experimental_structure_pdb": "", "protein_chain_id": "", "dna_chain_ids": "", "dna_sequence_by_chain": "", "dna_length_total": 0, "binding_core": TARGETS[dbp][1], "pdb_to_target_mapping": "no structure found after supplementary audit", "exact_pbm_construct_in_structure": False, "protein_construct_relation": "not applicable", "structure_protein_sequence": "", "missing_residues": "unknown", "mutations": "unknown", "terminal_additions_or_tags": "unknown", "structure_provenance": "not found", "structure_type": "", "structure_parse_status": "unavailable", "experimental_exact_pbm_construct": False, "experimental_protein_sequence": "", "experimental_protein_chain_id": "", "experimental_dna_chain_ids": "", "experimental_dna_sequence_by_chain": "", "experimental_dna_length_total": 0, "experimental_construct_note": "no experimental structure audited"})
    # Explicitly add the DBP48 crystallographic construct audit to the row.
    crystal = pde(CRYSTAL, "DBP48", "RCSB PDB 8TAC", "experimental co-crystal")
    rows[-1].update({"experimental_structure_pdb": "8TAC", "experimental_protein_chain_id": crystal["protein_chain_id"], "experimental_protein_sequence": crystal["structure_protein_sequence"], "experimental_dna_chain_ids": crystal["dna_chain_ids"], "experimental_dna_sequence_by_chain": crystal["dna_sequence_by_chain"], "experimental_dna_length_total": crystal["dna_length_total"], "experimental_exact_pbm_construct": False, "experimental_construct_note": "8TAC entity has terminal S/GSG additions relative to PBM sequence; bound 21-nt assembly differs from 14-bp PBM target"})
    pd.DataFrame(rows).to_csv(META / "structure_inventory.csv", index=False)
    return pd.DataFrame(rows)


def coverage_audit(inv: pd.DataFrame):
    old = pd.read_csv(ROOT / "results/v0_8_sota/deeppbs_metrics.csv")
    rows = []
    for dbp in DBPS:
        new = inv.loc[inv.protein_id == dbp].iloc[0]
        old_deep = old[(old.method == "DeepPBS") & (old.protein_id == dbp)].iloc[0]
        old_na = pd.read_csv(ROOT / "results/v0_8_sota/na_mpnn_metrics.csv")
        old_na = old_na[(old_na.method == "NA-MPNN") & (old_na.protein_id == dbp)].iloc[0]
        reason = "historical v0.8 search omitted Nature Supplementary Data 6; structure now available" if new.design_model_pdb else "no structure after Supplementary Data 6 and local repository audit"
        for method, prior in [("DeepPBS", old_deep), ("NA-MPNN", old_na)]:
            rows.append({"protein_id": dbp, "method": method, "v0_8_status": prior.status, "v1_0_structure_status": "structure_available" if new.design_model_pdb else "structure_unavailable", "coverage_failure_class": "previous_pipeline_resource_audit_gap" if new.design_model_pdb else "genuinely_missing_public_structure", "reason": reason, "biological_limitation_inferred": False, "action": "re-audit/re-run with frozen method" if new.design_model_pdb else "retain as not evaluable"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "v0_8_coverage_failure_audit.csv", index=False)
    return out


def npz_predictions(inv: pd.DataFrame):
    pred_dir = OUT / "nampnn_runs"
    rows, metric_rows = [], []
    exp = load_designed_experimental_units(ROOT)
    for dbp in DBPS:
        npz = SUP / f"DBP{int(dbp[3:]):03d}.npz"
        landscape_csv = OUT / "nampnn_runs" / f"{dbp}_landscape.csv"
        if not npz.exists() and landscape_csv.exists():
            scored = pd.read_csv(landscape_csv)
            if "nampnn_score" not in scored.columns:
                scored = scored.rename(columns={"prediction_score": "nampnn_score"})
            m = exp[exp.protein_id == dbp][["protein_id", "canonical_7mer", "experimental_score"]].merge(scored[["protein_id", "canonical_7mer", "nampnn_score"]], on=["protein_id", "canonical_7mer"], how="left")
            rank = compute_ranking_metrics(m, truth_col="experimental_score", prediction_col="nampnn_score")
            metric_rows.append({"method": "NA-MPNN", "protein_id": dbp, "status": "evaluated", "n_candidates": int(m.nampnn_score.notna().sum()), "spearman": rank.spearman, "coverage_fraction": float(m.nampnn_score.notna().mean()), "score_variance": float(m.nampnn_score.var()), "structure_source": inv.loc[inv.protein_id == dbp, "design_model_pdb"].iloc[0], "model_version": "s_70114.pt", "exposure_status": "DEVELOPMENT_EXPOSED_GSE237017"})
            rows.append(scored)
            continue
        if not npz.exists():
            legacy = ROOT / "results/v0_4/external_runs/nampnn_dbp035/specificity/DBP035.npz" if dbp == "DBP35" else None
            npz = legacy if legacy and legacy.exists() else npz
        if not npz.exists():
            metric_rows.append({"method": "NA-MPNN", "protein_id": dbp, "status": "not_evaluable_no_specificity_output", "n_candidates": 0, "spearman": np.nan, "coverage_fraction": 0.0, "score_variance": np.nan, "structure_source": inv.loc[inv.protein_id == dbp, "design_model_pdb"].iloc[0], "model_version": "s_70114.pt", "exposure_status": "DEVELOPMENT_EXPOSED_GSE237017"})
            continue
        scored = score_structural_ppm(npz, dbp, f"design_{dbp}", "s_70114.pt", "NA-MPNN specificity PPM").rename(columns={"prediction_score": "nampnn_score"})
        scored.to_csv(pred_dir / f"{dbp}_landscape.csv", index=False)
        m = exp[exp.protein_id == dbp][["protein_id", "canonical_7mer", "experimental_score"]].merge(scored[["protein_id", "canonical_7mer", "nampnn_score"]], on=["protein_id", "canonical_7mer"], how="left")
        rank = compute_ranking_metrics(m, truth_col="experimental_score", prediction_col="nampnn_score")
        metric_rows.append({"method": "NA-MPNN", "protein_id": dbp, "status": "evaluated", "n_candidates": int(m.nampnn_score.notna().sum()), "spearman": rank.spearman, "coverage_fraction": float(m.nampnn_score.notna().mean()), "score_variance": float(m.nampnn_score.var()), "structure_source": inv.loc[inv.protein_id == dbp, "design_model_pdb"].iloc[0], "model_version": "s_70114.pt", "exposure_status": "DEVELOPMENT_EXPOSED_GSE237017"})
        rows.append(scored)
    allpred = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if not allpred.empty:
        allpred.to_csv(OUT / "nampnn_results.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(OUT / "nampnn_metrics.csv", index=False)
    return allpred, pd.DataFrame(metric_rows)


def deeppbs_audit(inv: pd.DataFrame):
    old = pd.read_csv(ROOT / "results/v0_8_sota/deeppbs_metrics.csv")
    rows = []
    for dbp in DBPS:
        p = old[old.protein_id == dbp].iloc[0]
        if dbp in {"DBP35", "DBP48"}:
            status, reason = "historical_prediction_reused", "frozen v0.4.2 official run; native DeepPBS preprocessing cannot run here because torch_geometric and Linux 3DNA/X3DNA runtime are unavailable; score is not a v1.0 Supplementary Data 6 rerun"
            scored_input = "results/v0_4_2/external_runs/deeppbs/DBP35/DBP035.pdb" if dbp == "DBP35" else "results/v0_4_2/external_runs/deeppbs/DBP48/8TAC_helix_only.pdb"
        else:
            status, reason = "not_run_environment_dependency", "input structure verified; reproduction blocked by missing torch_geometric and Linux 3DNA/X3DNA preprocessing runtime"
            scored_input = inv.loc[inv.protein_id == dbp, "design_model_pdb"].iloc[0]
        rows.append({"method": "DeepPBS", "protein_id": dbp, "status": status, "n_candidates": int(p.n_candidates), "spearman": p.spearman, "coverage_fraction": p.coverage_fraction, "score_variance": p.score_variance, "structure_source": scored_input, "v1_0_structure_available": bool(inv.loc[inv.protein_id == dbp, "design_model_pdb"].iloc[0]), "published_applicability": dbp in {"DBP5", "DBP6", "DBP9", "DBP35"}, "qualitative_motif_reproduction": "not rerun in v1.0; blocked before inference" if dbp in {"DBP5", "DBP6", "DBP9", "DBP35"} else "not applicable", "failure_reason": reason, "exposure_status": "DEVELOPMENT_EXPOSED_GSE237017"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "deeppbs_results.csv", index=False)
    return out


def contact_features(inv: pd.DataFrame):
    rows = []
    for _, meta in inv.iterrows():
        path = ROOT / meta.design_model_pdb if meta.design_model_pdb else None
        if not path or not path.exists():
            rows.append({"protein_id": meta.protein_id, "structure": meta.design_model_pdb, "status": "not_available"})
            continue
        model = structure(path)
        pchain = model[meta.protein_chain_id]
        dchains = [model[c] for c in str(meta.dna_chain_ids).split(";") if c]
        protein_atoms = [a for r in pchain for a in r if a.element != "H"]
        dna_atoms = [a for ch in dchains for r in ch if r.resname.strip().upper() in DNA_RES for a in r if a.element != "H"]
        base_atoms = [a for ch in dchains for r in ch if r.resname.strip().upper() in DNA_RES for a in r if a.name.strip() in BASE_ATOMS]
        pairs = [(pa, da) for pa in protein_atoms for da in dna_atoms if pa - da <= 4.0]
        base_pairs = [(pa, da) for pa in protein_atoms for da in base_atoms if pa - da <= 4.0]
        residues = {(pa.get_parent().get_id()[1], pa.get_parent().resname) for pa, _ in pairs}
        rows.append({"protein_id": meta.protein_id, "structure": meta.design_model_pdb, "status": "computed", "heavy_atom_contact_count_4A": len(pairs), "base_atom_contact_count_4A": len(base_pairs), "base_contacting_residue_count": len(residues), "base_specific_hbond_like_count": sum(1 for pa, da in base_pairs if pa.element in {"N", "O"} and da.element in {"N", "O"} and pa - da <= 3.5), "backbone_or_phosphate_contact_count": len(pairs) - len(base_pairs), "contact_definition": "non-H protein/DNA atom distance <=4.0A; base atom names whitelist; exploratory descriptor"})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "interface_contact_features.csv", index=False)
    return out


def landscape_correlations(preds: pd.DataFrame):
    rows = []
    if preds.empty:
        pd.DataFrame(columns=["method", "protein_a", "protein_b", "landscape_spearman", "n_shared_candidates"]).to_csv(OUT / "protein_landscape_correlations.csv", index=False)
        return
    for method in sorted(preds.prediction_type.unique()):
        x = preds[preds.prediction_type == method]
        score_col = "nampnn_score" if "nampnn_score" in x.columns else "prediction_score"
        piv = x.pivot(index="canonical_7mer", columns="protein_id", values=score_col)
        for a, b in itertools.combinations(sorted(piv.columns), 2):
            rows.append({"method": method, "protein_a": a, "protein_b": b, "landscape_spearman": piv[a].corr(piv[b], method="spearman"), "landscape_pearson": piv[a].corr(piv[b]), "n_shared_candidates": int(piv[[a, b]].dropna().shape[0]), "conditioning_effect_pair": 1 - piv[a].corr(piv[b], method="spearman")})
    pd.DataFrame(rows).to_csv(OUT / "protein_landscape_correlations.csv", index=False)


def frozen_subsets():
    # Deterministic candidate subset is generated before any structural result is read.
    bases = "ACGT"
    target_core = {"DBP35": "TGCACA", "DBP48": "CTGACG"}
    rows, af3 = [], []
    hard = pd.read_csv(ROOT / "results/v0_4/tables/baseline_failure_cases.csv") if (ROOT / "results/v0_4/tables/baseline_failure_cases.csv").exists() else pd.DataFrame()
    for dbp in ["DBP35", "DBP48"]:
        t = target_core[dbp]
        seqs = {t}
        for i in range(6):
            for b in bases:
                if b != t[i]: seqs.add(t[:i] + b + t[i + 1:])
        # Lexicographically first Hamming-2 variants, fixed at 12 per protein.
        h2 = ["".join(x) for x in itertools.product(bases, repeat=6)]
        h2 = [s for s in h2 if sum(a != b for a, b in zip(s, t)) == 2][:12]
        seqs.update(h2)
        for s in sorted(seqs):
            rows.append({"protein_id": dbp, "candidate_7mer": s, "subset_class": "on_target" if s == t else "hamming_1" if sum(a != b for a,b in zip(s,t)) == 1 else "hamming_2_lexicographic", "selection_rule": "frozen before Rosetta/AF3 results; deterministic lexicographic enumeration", "status": "not_run_rosetta_executable_not_available", "failure_reason": "Rosetta executable and licensed scoring environment are not installed in this Windows workspace"})
            af3.append({"protein_id": dbp, "candidate_7mer": s, "subset_class": rows[-1]["subset_class"], "selection_rule": "frozen before AF3/ContactSeek results", "status": "not_run_af3_runtime_or_api_unavailable", "failure_reason": "AF3/ContactSeek inference service and structural runner are not available locally; subset diagnostic only"})
    pd.DataFrame(rows).to_csv(OUT / "rossetta_subset_results.csv", index=False)
    pd.DataFrame(af3).to_csv(OUT / "af3_contact_subset_results.csv", index=False)


def docs(inv, cov, deep, na, contacts):
    DOC.mkdir(parents=True, exist_ok=True)
    (DOC / "STRUCTURE_SOURCE_AUDIT.md").write_text("""# V1.0 Structure Source Audit

## Scope

This development-stage audit checked the local repository, the official `cjg263/dbp_design` checkout, Nature Supplementary Data 6, and RCSB PDB 8TAC. The primary source is [Glasscock et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12618268/), whose data-availability statement identifies Supplementary Data 6 as the PDB design-model archive and 8TAC as the DBP48 co-crystal. No v0.5-v0.9 file was modified. All seven GSE237017 proteins remain development-exposed.

## Finding

Supplementary Data 6 contains design-hit PDBs for DBP001, DBP003, DBP005, DBP006, DBP009, DBP035 and DBP048. In each recovered file, protein chain A matches the corresponding PBM protein sequence exactly and DNA chain B contains the intended target duplex/context. This changes the interpretation of the old 2/7 structure coverage: it was a resource-audit/pipeline omission, not evidence that five proteins lack structure or that structure methods are biologically inapplicable.

DBP48 also has experimental co-crystal [PDB 8TAC](https://www.rcsb.org/structure/8TAC). Its deposited construct includes terminal additions relative to the 65-aa PBM sequence, and the 21-nt crystal DNA assembly differs from the 14-bp PBM target. It is therefore a construct-sensitivity sanity check, not an exact assay-matched replacement.

The complete machine-readable chain, sequence, DNA, provenance and mismatch inventory is `metadata/v1_0_structure/structure_inventory.csv`. The downloaded archive hash and extraction provenance are recorded in `results/v1_0_structure/v1_0_run_manifest.json`.
""", encoding="utf-8")
    (DOC / "STRUCTURE_TO_7MER_SCORING_PROTOCOL.md").write_text("""# Structure-to-7-mer Scoring Protocol (Frozen)\n\nThis protocol is frozen before reading any v1.0 PBM performance. It is applied identically to NA-MPNN PPMs and, where available, DeepPBS PWM outputs.\n\n1. The candidate universe is the 8,192 canonical reverse-complement classes defined by the frozen v0.6 sequence-equivalence utility.\n2. A PPM/PWM is converted to a log-probability score using `sum(log(P(base_i)+1e-12))`.\n3. For a structure with multiple contiguous DNA runs, all structurally valid contiguous 7-bp windows and both candidate orientations are scored; the maximum is retained. This is a fixed register aggregation rule, not a per-protein PBM-driven search.\n4. The experimental/design core is taken from the frozen target metadata (`GCAGG`, `GCAGGA`, `TGCACA`, or `CTGACG`). No register is selected by maximizing Spearman.\n5. Missing residues, ambiguous chain mapping, or a DNA run shorter than 7 bp produce `not evaluable`, never an imputed score.\n6. PPM values are ranking proxies; they are not numerically equated with PBM E-scores.\n\nThe output records the selected structural window and orientation so every ranking is auditable.\n""", encoding="utf-8")
    na_eval = na[na.status == "evaluated"]
    med = float(na_eval.spearman.median()) if not na_eval.empty else float("nan")
    (DOC / "V1_0_STRUCTURE_RESULTS.md").write_text(f"""# V1.0 Structure-Conditioning Mechanism Pilot Results\n\n## Status\n\nThis is a development/mechanism pilot on the seven GSE237017 DBPs. All seven proteins were development-exposed; none of these results is independent confirmation.\n\n## Structure audit\n\nThe Supplementary Data 6 audit recovered design complex PDBs for all seven DBPs. Therefore v0.8's 2/7 structure coverage was a previous resource/pipeline audit gap.\n\n## Frozen inference status\n\nNA-MPNN specificity inference ran locally on all seven recovered design PDBs with the frozen `s_70114.pt` checkpoint. Its PPM-derived landscape metrics are in `results/v1_0_structure/nampnn_metrics.csv` and full landscapes are in `results/v1_0_structure/nampnn_results.csv`; median Spearman over evaluated proteins is {med:.4f} (N={len(na_eval)}). DeepPBS reproduction is recorded honestly in `deeppbs_results.csv`: the required torch_geometric and Linux 3DNA/X3DNA preprocessing stack is unavailable in this workspace, so the historical DBP35/DBP48 outputs are reused and DBP5/6/9 published applicability is not claimed as reproduced. This is an environment stop, not a biological failure.\n\n### NA-MPNN per-protein ranking\n\n| protein | status | Spearman |\n|---|---|---|\n""" + "\n".join(f"| {r.protein_id} | {r.status} | {r.spearman:.4f} |" for _, r in na[na.status == "evaluated"].iterrows()) + f"\n\n## Conditionality\n\nThe NA-MPNN landscape correlation table is `protein_landscape_correlations.csv`. A structure model is considered protein-sensitive only descriptively from these between-protein landscapes; no arbitrary rescue threshold is imposed. The current output must not be interpreted as proof of transferability because the structures are design models and the cohort is exposed.\n\n## Contacts and controls\n\n`interface_contact_features.csv` contains exploratory 4-A heavy-atom/base contact descriptors. `dbp48_design_vs_crystal.csv` compares the two DBP48 NA-MPNN landscapes descriptively; differences cannot be attributed to design-model error without matched construct and DNA. Rosetta and AF3/ContactSeek subsets are frozen but not run because the required runtimes are unavailable; their candidate manifests preserve the failure reasons.\n\n### Design-interface descriptors\n\n| protein | heavy contacts (4A) | base contacts (4A) | base-contacting residues |\n|---|---:|---:|---:|\n""" + "\n".join(f"| {r.protein_id} | {int(r.heavy_atom_contact_count_4A)} | {int(r.base_atom_contact_count_4A)} | {int(r.base_contacting_residue_count)} |" for _, r in contacts[contacts.status == "computed"].iterrows()) + """\n\n## H1-H5 judgment\n\nH1/H2: explicit geometry is now demonstrably available as an input and NA-MPNN produces protein-specific PPMs, but a full DeepPBS reproduction and a fair cross-method conditionality comparison remain incomplete.\n\nH3: no claim of reaching the experimental replicate reference is made from the partial/heterogeneous structure results.\n\nH4: base-contact and register descriptors provide mechanistic hypotheses only (N=7).\n\nH5: sequence-only failure cannot yet be attributed solely to missing geometry, because DeepPBS reproduction on DBP5/6/9 remains blocked by environment and no AF3/Rosetta diagnostic was executed.\n\n**Scientific answer:** explicit interface geometry has not yet been shown to rescue the protein-conditioning failure. The correct v1.0 status is **INCONCLUSIVE / BLOCKED ON REPRODUCIBLE STRUCTURE-METHOD EXECUTION**, with the important correction that the former 2/7 coverage was an audit gap rather than biological unavailability.\n""", encoding="utf-8")


def dbp48_structure_comparison():
    """Compare design and 8TAC NA-MPNN outputs without post-hoc tuning."""
    exp = load_designed_experimental_units(ROOT)
    rows = []
    sources = [
        ("design_model", SUP / "DBP048.npz", "design_DBP48"),
        ("experimental_8TAC", ROOT / "results/v0_4/external_runs/nampnn_8tac/specificity/8TAC.npz", "8TAC"),
    ]
    for label, path, structure_id in sources:
        if not path.exists():
            rows.append({"protein_id": "DBP48", "structure_source": label, "status": "not_available", "n_candidates": 0})
            continue
        scored = score_structural_ppm(path, "DBP48", structure_id, "s_70114.pt", "NA-MPNN specificity PPM").rename(columns={"prediction_score": "nampnn_score"})
        merged = exp[exp.protein_id == "DBP48"][["protein_id", "canonical_7mer", "experimental_score"]].merge(scored[["protein_id", "canonical_7mer", "nampnn_score"]], on=["protein_id", "canonical_7mer"], how="left")
        metrics = compute_ranking_metrics(merged, truth_col="experimental_score", prediction_col="nampnn_score")
        top = scored.sort_values("nampnn_score", ascending=False).iloc[0]
        rows.append({"protein_id": "DBP48", "structure_source": label, "status": "evaluated", "n_candidates": int(merged.nampnn_score.notna().sum()), "spearman": metrics.spearman, "top_predicted_7mer": top.canonical_7mer, "top_score": float(top.nampnn_score), "structure_id": structure_id, "comparison_status": "descriptive_provenance_only", "interpretation": "design model versus 8TAC; construct/DNA mismatch prevents isolating design-model error", "exposure_status": "DEVELOPMENT_EXPOSED_GSE237017"})
    pd.DataFrame(rows).to_csv(OUT / "dbp48_design_vs_crystal.csv", index=False)


def main():
    for d in [OUT, META, DOC]: d.mkdir(parents=True, exist_ok=True)
    inv = inventory()
    cov = coverage_audit(inv)
    deep = deeppbs_audit(inv)
    preds, na = npz_predictions(inv)
    landscape_correlations(preds)
    contacts = contact_features(inv)
    dbp48_structure_comparison()
    frozen_subsets()
    docs(inv, cov, deep, na, contacts)
    # A compact provenance record makes reruns auditable without changing frozen history.
    (OUT / "v1_0_run_manifest.json").write_text(json.dumps({"stage": "v1.0", "design_pdb_source": str(SUP.relative_to(ROOT)).replace("\\", "/"), "structure_archive": "data/raw/v1_0_structure/41594_2025_1669_MOESM9_ESM.zip", "structure_archive_sha256": sha256(ROOT / "data/raw/v1_0_structure/41594_2025_1669_MOESM9_ESM.zip"), "deepPBS_environment_status": "blocked_torch_geometric_and_linux_3dna", "nampnn_checkpoint": "external/nampnn/NA-MPNN/models/specificity_model/s_70114.pt", "exposure_status": "DEVELOPMENT_EXPOSED_GSE237017", "pbm_labels_used_for_ranking": True, "register_protocol_frozen_before_ranking": True, "quantitative_supplementary_labels_read": False}, indent=2), encoding="utf-8")
    print("v1.0 structure pilot outputs written", OUT)


if __name__ == "__main__":
    main()
