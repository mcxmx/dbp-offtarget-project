from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from src.utils import ensure_dir, project_root


ROOT = project_root()
DOCS = ensure_dir(ROOT / "docs" / "v0_4")
METADATA = ensure_dir(ROOT / "metadata" / "v0_4")
RESULTS = ensure_dir(ROOT / "results" / "v0_4")
TABLES = ensure_dir(RESULTS / "tables")
LOGS = ensure_dir(ROOT / "logs" / "v0_4")


def git_head(path: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
    except Exception as exc:  # pragma: no cover - provenance fallback
        return f"unavailable: {exc}"


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(col) for col in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in frame.columns) + " |")
    return "\n".join(lines)


def discover_dbp35_structure() -> dict:
    pdb_path = ROOT / "external" / "dbp_design" / "2b_design_mpnn" / "DBP035.pdb"
    dna_seq = "unknown"
    protein_len = None
    dna_len = None
    if pdb_path.exists():
        protein_residues = []
        dna_residues = []
        last_key = None
        residue_to_base = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "A": "A", "C": "C", "G": "G", "T": "T"}
        for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21].strip()
            resnum = line[22:26].strip()
            resname = line[17:20].strip()
            key = (chain, resnum, resname)
            if key == last_key:
                continue
            last_key = key
            if chain == "A":
                protein_residues.append(key)
            elif chain == "B":
                dna_residues.append(key)
        dna_seq = "".join(residue_to_base.get(resname, "N") for _, _, resname in dna_residues)
        protein_len = len(protein_residues)
        dna_len = len(dna_residues)
    return {
        "protein_id": "DBP35",
        "structure_type": "designed complex model",
        "experimental_or_predicted": "predicted/theoretical Rosetta model",
        "pdb_id": pd.NA,
        "local_file": "external/dbp_design/2b_design_mpnn/DBP035.pdb" if pdb_path.exists() else pd.NA,
        "source_url": "https://github.com/cjg263/dbp_design",
        "paper_reference": "Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4",
        "bound_dna_sequence": dna_seq,
        "bound_dna_length": dna_len,
        "protein_chain": "A",
        "dna_chain": "B",
        "protein_length_in_structure": protein_len,
        "structure_confidence": "medium",
        "notes": "Only readily available designed-DBP complex model found in the official dbp_design repository checkout. Header marks EXPDTA THEORETICAL MODEL / Rosetta. DNA chain B contains target-B duplex context; NA-MPNN featurizer masks 29 of 30 DNA residues.",
    }


def discover_dbp48_structure() -> dict:
    cif_path = ROOT / "data" / "raw" / "rcsb" / "mmcif" / "8TAC.cif"
    return {
        "protein_id": "DBP48",
        "structure_type": "experimental protein-DNA complex",
        "experimental_or_predicted": "experimental X-ray diffraction",
        "pdb_id": "8TAC",
        "local_file": "data/raw/rcsb/mmcif/8TAC.cif" if cif_path.exists() else pd.NA,
        "source_url": "https://www.rcsb.org/structure/8TAC",
        "paper_reference": "Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4; PMID 40940539",
        "bound_dna_sequence": "ACCTGACGCGA;TTCGCGTCAGG",
        "bound_dna_length": 21,
        "protein_chain": "B",
        "dna_chain": "C;D",
        "protein_length_in_structure": 66,
        "structure_confidence": "high",
        "notes": "RCSB mmCIF reports X-RAY DIFFRACTION at 2.34 A and chains A/B as the designed protein entity; v0.2 curation selected protein chain B and DNA chain C. NA-MPNN is run on the full biological assembly, so this is a structure-level diagnostic, not isolated single-chain scoring.",
    }


def make_structure_manifest() -> pd.DataFrame:
    protein_ids = ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]
    dbp35 = discover_dbp35_structure()
    dbp48 = discover_dbp48_structure()
    rows = []
    for protein_id in protein_ids:
        if protein_id == "DBP35":
            rows.append(dbp35)
        elif protein_id == "DBP48":
            rows.append(dbp48)
        else:
            rows.append(
                {
                    "protein_id": protein_id,
                    "structure_type": "not found in current official checkout",
                    "experimental_or_predicted": "unknown/not available",
                    "pdb_id": pd.NA,
                    "local_file": pd.NA,
                    "source_url": "https://github.com/cjg263/dbp_design",
                    "paper_reference": "Computational design of sequence-specific DNA-binding proteins; DOI 10.1038/s41594-025-01669-4",
                    "bound_dna_sequence": pd.NA,
                    "bound_dna_length": pd.NA,
                    "protein_chain": pd.NA,
                    "dna_chain": pd.NA,
                    "protein_length_in_structure": pd.NA,
                    "structure_confidence": "none",
                    "notes": "No protein-DNA complex structure/model for this designed DBP was found in the checked official dbp_design repository. Structure-aware baselines are not evaluated for this protein in v0.4.",
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(METADATA / "designed_dbp_structure_manifest.csv", index=False)
    return df


def _split_hits() -> dict[tuple[str, str], list[str]]:
    hits: dict[tuple[str, str], list[str]] = {}
    split_sources = {
        "DeepPBS": list((ROOT / "external" / "deeppbs" / "DeepPBS" / "run" / "folds").glob("*.txt")),
        "NA-MPNN": list((ROOT / "external" / "nampnn" / "NA-MPNN" / "splits").glob("*.json")),
    }
    tokens = ["8tac", "dbp1", "dbp3", "dbp5", "dbp6", "dbp9", "dbp35", "dbp48", "gse237017"]
    for baseline, paths in split_sources.items():
        for path in paths:
            text = path.read_text(encoding="utf-8", errors="replace").lower()
            for token in tokens:
                if token in text:
                    hits.setdefault((baseline, token), []).append(path.name)
    return hits


def make_overlap_audit() -> pd.DataFrame:
    protein_ids = ["DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"]
    baselines = ["DeepPBS", "NA-MPNN"]
    rows = []
    hits = _split_hits()
    motifs = pd.read_csv(ROOT / "metadata" / "v0_3_1" / "designed_dbp_target_definitions.csv")
    for protein_id in protein_ids:
        motif = motifs.loc[motifs["protein_id"] == protein_id, "designed_binding_site_motif"].iloc[0]
        for baseline in baselines:
            exact = False
            structure_seen = False
            motif_seen = False
            risk_level = "low_to_medium_unresolved_homology"
            evidence = "Direct scan of checked-out official split files found no DBP/GSE direct hit for this protein. Homolog-level leakage is not excluded."
            if baseline == "NA-MPNN" and protein_id == "DBP48" and ("NA-MPNN", "8tac") in hits:
                exact = True
                structure_seen = True
                risk_level = "high_not_zero_shot"
                evidence = (
                    "8tac appears in NA-MPNN split files: "
                    + ";".join(hits[("NA-MPNN", "8tac")])
                    + ". 8TAC is the DBP48 experimental structure used here."
                )
            rows.append(
                {
                    "protein_id": protein_id,
                    "baseline": baseline,
                    "exact_protein_seen": bool(exact),
                    "homolog_seen": "not_assessed_by_sequence_search",
                    "structure_seen": bool(structure_seen),
                    "motif_or_target_seen": bool(motif_seen),
                    "evidence": evidence,
                    "risk_level": risk_level,
                    "notes": "A full homolog-level leakage audit requires model training-set protein sequences, which are not distributed in the checked repositories. Do not call any result strict zero-shot until this is resolved.",
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(METADATA / "baseline_data_overlap_audit.csv", index=False)
    return df


def make_external_baseline_provenance() -> pd.DataFrame:
    rows = [
        {
            "baseline": "DeepPBS",
            "repository_url": "https://github.com/timkartar/DeepPBS",
            "local_path": "external/deeppbs/DeepPBS",
            "commit_hash": git_head(ROOT / "external" / "deeppbs" / "DeepPBS"),
            "weights_or_model_files": "external/deeppbs/DeepPBS/run/output/*/Model.best.tar; bundled ensemble listed under run/plot_scripts/txts/",
            "inference_status_v0_4": "not_evaluable_in_current_environment",
            "status_reason": "Official preprocessing requires Linux-oriented 3DNA/Curves-style process chain plus torch_geometric/freesasa/pdb2pqr setup. Current environment is Windows PowerShell without conda/docker. No fair DeepPBS prediction is generated in v0.4.",
        },
        {
            "baseline": "NA-MPNN",
            "repository_url": "https://github.com/baker-laboratory/NA-MPNN",
            "local_path": "external/nampnn/NA-MPNN",
            "commit_hash": git_head(ROOT / "external" / "nampnn" / "NA-MPNN"),
            "weights_or_model_files": "external/nampnn/NA-MPNN/models/specificity_model/s_70114.pt",
            "inference_status_v0_4": "ran_for_dbp35_and_dbp48_only",
            "status_reason": "Official Windows CPU inference ran on the DBP035 theoretical design model and RCSB 8TAC. DBP48/8TAC appears in NA-MPNN validation split files and is not a zero-shot designed-DBP result. The other five designed DBPs lack public complex structures in the current checkout.",
        },
        {
            "baseline": "SimpleProteinConditionalBaseline",
            "repository_url": "local project scaffold",
            "local_path": "src/models/simple_protein_conditional_baseline.py",
            "commit_hash": "local",
            "weights_or_model_files": "none",
            "inference_status_v0_4": "not_trained_not_evaluated",
            "status_reason": "No assay-compatible natural PBM/uPBM training set has been incorporated. v0.4 does not train on designed DBP rows as a main generalization result.",
        },
    ]
    df = pd.DataFrame(rows)
    df.to_csv(METADATA / "external_baseline_provenance.csv", index=False)
    return df


def write_docs(structure_manifest: pd.DataFrame) -> None:
    feasibility = "\n".join(
        [
            "# v0.4 Baseline Feasibility Audit",
            "",
            "Audit date: 2026-09-02",
            "",
            "## Benchmark Held Fixed",
            "",
            "Input benchmark is v0.3.1: 7 designed DBPs, 57,344 protein-RC-class experimental units, primary score `experimental_escore_consensus` from processed GSE237017 uPBM E-scores. v0.3.1 files are not modified.",
            "",
            "## DeepPBS",
            "",
            "- Official code: https://github.com/timkartar/DeepPBS",
            "- Local path: `external/deeppbs/DeepPBS`",
            f"- Commit: `{git_head(ROOT / 'external' / 'deeppbs' / 'DeepPBS')}`",
            "- Paper DOI from README: 10.1038/s41592-024-02372-w",
            "- Official output: predicted DNA position weight matrix / PWM from a protein-DNA complex structure.",
            "- Required input: biological protein-DNA assembly as PDB/CIF, then DeepPBS preprocessing into graph/shape features.",
            "- Structure requirement: yes. The method is not a raw `(protein sequence, DNA 7-mer)` scorer.",
            "- Fixed cognate structure to alternative DNA ranking: only possible after a justified PWM-to-7mer mapping; it does not directly predict PBM E-score.",
            "- Designed DBP overlap: no direct GSE237017/DBP ID or 8TAC hit in checked DeepPBS fold files; full homolog overlap cannot be excluded without the complete training-set sequence list.",
            "- v0.4 execution status: not evaluated. The official preprocessing chain is Linux-oriented and requires 3DNA/Curves-style processing plus `torch_geometric`, `freesasa`, and related dependencies. Current Windows environment is not a fair execution target.",
            "",
            "## NA-MPNN",
            "",
            "- Official code: https://github.com/baker-laboratory/NA-MPNN",
            "- Local path: `external/nampnn/NA-MPNN`",
            f"- Commit: `{git_head(ROOT / 'external' / 'nampnn' / 'NA-MPNN')}`",
            "- Official output: predicted PPM over residue types for nucleic-acid positions from a protein-nucleic-acid structure.",
            "- Required input: protein-DNA/RNA structure.",
            "- Structure requirement: yes. The specificity mode is not a structure-free `(protein sequence, DNA 7-mer)` scorer.",
            "- Fixed cognate structure to alternative DNA ranking: possible only by mapping PPM log-probabilities over a structure-defined DNA window into candidate 7-mer scores. This is a derived task-harmonization layer.",
            "- Designed DBP overlap: `8tac` appears in NA-MPNN `design_valid.json` and `specificity_valid.json`; DBP48/8TAC is therefore not zero-shot. No DBP35/GSE237017 direct hit was found in checked split files, but homolog-level risk is unresolved.",
            "- v0.4 execution status: official inference ran successfully on the bundled 1am9 example, on `external/dbp_design/2b_design_mpnn/DBP035.pdb`, and on RCSB `8TAC`.",
            "",
            "## TransBind",
            "",
            "TransBind is not executed in v0.4. It remains a task-comparison item until its input/output can be fairly mapped to this processed uPBM 7-mer ranking benchmark.",
            "",
            "## Structure Availability Summary",
            "",
            markdown_table(structure_manifest),
            "",
            "## Main Feasibility Conclusion",
            "",
            "v0.4 can fairly evaluate sequence-only baselines on all seven proteins and partial NA-MPNN structural-PPM baselines on DBP35 and DBP48. DeepPBS is recorded as not evaluable in this environment. DBP48/8TAC has explicit NA-MPNN split overlap risk and is not zero-shot.",
            "",
        ]
    )
    write_text(DOCS / "BASELINE_FEASIBILITY_AUDIT.md", feasibility)
    write_text(
        DOCS / "TASK_HARMONIZATION.md",
        """# v0.4 Task Harmonization\n\nThe ground truth task is per-protein ranking of reverse-complement canonical 7-mers by processed experimental uPBM E-score consensus.\n\n## PBM E-score\n\n`experimental_escore_consensus` is a processed experimental uPBM E-score consensus. It is suitable here as a per-protein ranking target. It is not Kd, binding free energy, binding probability, or an absolute cross-protein affinity scale.\n\n## Sequence-Only Baselines\n\nHamming/edit/k-mer metrics are sequence-only proxy metrics computed against the v0.3.1 paper motif reference. They are not protein-conditioned.\n\n## DeepPBS\n\nDeepPBS predicts a structure-conditioned DNA PWM. A fair mapping to this benchmark would require a validated protein-DNA complex structure for each designed DBP and a predeclared PWM-to-7mer scoring rule. v0.4 does not generate DeepPBS predictions because the official preprocessing chain is not runnable in the current environment.\n\n## NA-MPNN\n\nNA-MPNN specificity mode predicts a PPM over residue types at nucleic-acid positions in a supplied structure. For DBP35 and DBP48, v0.4 maps the predicted DNA-position probabilities to canonical 7-mer scores as follows:\n\n1. Use the official `s_70114.pt` specificity checkpoint.\n2. Extract predicted probabilities for DA/DC/DG/DT at nucleic-acid positions from the official inference `.npz` output.\n3. Preserve each contiguous DNA-chain run rather than concatenating unrelated chains.\n4. For each candidate 7-mer, compute the sum of log probabilities in every complete 7-position window.\n5. Score a canonical RC class by the maximum score over candidate orientation and reverse-complement orientation.\n\nThe resulting score is named `partial_structural_ppm_best_window_log_probability`. It is a derived ranking score, not a PBM E-score and not an affinity. DBP48/8TAC has NA-MPNN validation-split overlap and must be interpreted as a diagnostic, not a zero-shot result.\n\n## Main Metrics\n\nMetrics are computed per protein and macro-summarized. The primary ranking metric is Spearman correlation; NDCG@1%, NDCG@5%, top-1% recovery, and sampled pairwise ranking accuracy are secondary. Rows from different proteins are not pooled into a single main correlation.\n""",
    )
    write_text(
        DOCS / "FAILURE_RESOLUTION_DEFINITION.md",
        """# v0.4 Failure Resolution Definition\n\nThis document is predeclared before inspecting v0.4 failure counts.\n\n## Starting Set\n\nUse v0.3.1 sequence-vs-experiment disagreement candidates: per protein, processed uPBM E-score consensus in the top 5% and RC-aware Hamming similarity to the paper motif at or below the per-protein median.\n\n## Resolution Rule\n\nFor any protein-conditioned baseline with predictions on that protein, a v0.3.1 disagreement candidate is considered resolved when the baseline ranks it in the top 10% of that protein's predicted scores.\n\nThis asks whether a protein-conditioned method elevates experimentally high-scoring sequences that a sequence-only Hamming proxy did not prioritize. It does not imply the method correctly models physical binding mechanism.\n\n## Not-Evaluable Cases\n\nIf a baseline has no prediction for a protein, candidates from that protein are counted as not evaluable for that baseline, not unresolved.\n""",
    )
    write_text(
        DOCS / "NATURAL_PBM_CONTROL_PLAN.md",
        """# Natural PBM Control Plan for v0.5+\n\nThe v0.4 designed-DBP benchmark uses processed uPBM E-scores from GSE237017. A future natural-to-designed OOD claim needs assay-matched controls because a natural HT-SELEX to designed uPBM comparison would confound protein distribution shift with assay distribution shift.\n\nPriority sources to investigate next:\n\n- CIS-BP PBM-derived specificity tables where protein IDs and protein sequences are resolvable.\n- UniPROBE / universal PBM datasets with raw or processed k-mer enrichment and protein identity.\n- JASPAR entries backed by PBM/uPBM experiments, used only when the experimental provenance is clear.\n\nRecommended comparisons:\n\n- Natural PBM train to natural PBM held-out protein.\n- Natural PBM train to designed GSE237017 uPBM external test.\n- Protein-family-out splits within natural PBM.\n- Designed leave-one-protein/cluster-out only as secondary, because N=7 is small.\n\nAll DNA splits must group reverse-complement equivalent sequences together.\n""",
    )
    write_text(
        DOCS / "PROPOSED_MODEL_REQUIREMENTS.md",
        """# Proposed Model Requirements Placeholder\n\nThis file is finalized by `analysis/v0_4/05_write_reports.py` after baseline and failure analyses are generated. v0.4 does not implement a new proposed model.\n""",
    )


def main() -> None:
    structure_manifest = make_structure_manifest()
    make_overlap_audit()
    make_external_baseline_provenance()
    write_docs(structure_manifest)
    pip_freeze = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
    (LOGS / "v0_4_python_environment_pip_freeze.txt").write_text(pip_freeze, encoding="utf-8")
    summary = {
        "deeppbs_commit": git_head(ROOT / "external" / "deeppbs" / "DeepPBS"),
        "nampnn_commit": git_head(ROOT / "external" / "nampnn" / "NA-MPNN"),
        "dbp_design_commit": git_head(ROOT / "external" / "dbp_design"),
        "n_structures_found": int(structure_manifest["local_file"].notna().sum()),
        "structures_found_for": structure_manifest.loc[structure_manifest["local_file"].notna(), "protein_id"].tolist(),
    }
    (METADATA / "v0_4_feasibility_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
