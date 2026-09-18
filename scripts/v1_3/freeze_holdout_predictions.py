from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from file_manifest import record_file


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/v1_3/holdout_predictions_unscored.tsv"
REPORT = ROOT / "reports/v1.3/HOLDOUT_PREDICTION_FREEZE.md"
BASES = "ACGT"
HOLDOUTS = ["DBP023", "DBP056", "DBP062"]
EPS = 1e-12
PDB_DIR = ROOT / "data/raw/v1_0_structure/design_pdbs/design_pdbs"


def pdb_dna(protein: str) -> tuple[str, str, str]:
    path = PDB_DIR / f"{protein}.pdb"
    residues, seen = [], set()
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        residue, chain = line[17:20].strip(), line[21].strip()
        key = (chain, line[22:27])
        if chain == "B" and residue in {"DA", "DC", "DG", "DT"} and key not in seen:
            seen.add(key)
            residues.append((int(line[22:26]), residue[-1]))
    sequence = "".join(base for _, base in residues)
    if len(sequence) % 2 or sequence[len(sequence) // 2:] != reverse_complement(sequence[: len(sequence) // 2]):
        raise ValueError(f"{protein}: DNA chain B is not a forward strand followed by its reverse complement")
    first = sequence[: len(sequence) // 2]
    return first, str(residues[0][0]), str(residues[len(first) - 1][0])


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze DeepPBS holdout predictions without opening competition labels.")
    parser.add_argument("--prediction-dir", type=Path, required=True)
    args = parser.parse_args()
    rows, summaries = [], []
    for protein in HOLDOUTS:
        path = args.prediction_dir / f"{protein}.npz_predict.npz"
        if not path.exists():
            raise FileNotFoundError(path)
        data = np.load(path, allow_pickle=False)
        P, seq_onehot = np.asarray(data["P"], float), np.asarray(data["Seq"], float)
        if P.ndim != 2 or P.shape[1] != 4 or seq_onehot.shape != P.shape:
            raise ValueError(f"{protein}: unexpected P/Seq shapes {P.shape}/{seq_onehot.shape}")
        if not np.isfinite(P).all() or (P <= 0).any() or not np.allclose(P.sum(axis=1), 1, atol=1e-6):
            raise ValueError(f"{protein}: invalid probability matrix")
        sequence = "".join(BASES[i] for i in np.argmax(seq_onehot, axis=1))
        pdb_sequence, dna_start, dna_end = pdb_dna(protein)
        if sequence != pdb_sequence:
            raise ValueError(f"{protein}: DeepPBS sequence {sequence} != frozen PDB first strand {pdb_sequence}")
        for position, wt in enumerate(sequence, start=1):
            for mutant in BASES:
                if mutant == wt:
                    continue
                effect = math.log(P[position - 1, BASES.index(mutant)] + EPS) - math.log(P[position - 1, BASES.index(wt)] + EPS)
                rows.append({
                    "protein": protein,
                    "position": position,
                    "wt_base": wt,
                    "mutant_base": mutant,
                    "predicted_effect": effect,
                    "method": "DeepPBS_native_pwm_frozen_ensemble",
                    "structure_source": f"data/raw/v1_0_structure/design_pdbs/design_pdbs/{protein}.pdb",
                    "protein_chain": "A",
                    "dna_chain": "B",
                    "dna_residue_span": f"{dna_start}-{dna_end}",
                    "structure_first_strand": sequence,
                    "mapping_rule": "position i maps to DeepPBS/PDB first-strand index i; no register search",
                    "orientation": "PDB chain B first strand as emitted by DeepPBS",
                    "score_definition": "log(P_mutant+1e-12)-log(P_wt+1e-12)",
                    "deeppbs_commit": "8bfb211dd67f02877841f6f33aa493ddf7daedf9",
                    "checkpoint_list": "DeepPBS.txt bundled five-checkpoint ensemble",
                    "prediction_source": path.as_posix(),
                })
        summaries.append({
            "protein": protein,
            "positions": len(sequence),
            "rows": 3 * len(sequence),
            "sequence": sequence,
            "prediction_min": float(np.min([r["predicted_effect"] for r in rows if r["protein"] == protein])),
            "prediction_max": float(np.max([r["predicted_effect"] for r in rows if r["protein"] == protein])),
        })
    out = pd.DataFrame(rows).sort_values(["protein", "position", "mutant_base"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    manifest = record_file(OUT)
    summary = pd.DataFrame(summaries)
    report = f"""# Holdout Prediction Freeze

Freeze status: **PREDICTIONS FROZEN; LABELS UNREAD**. Date: 2026-09-17.

This file records Stage 1 of the two-stage prospective evaluation. `scripts/v1_3/freeze_holdout_predictions.py` accepts only DeepPBS NPZ predictions and design PDBs. It does not import, open, or join any competition workbook or processed competition table.

## Frozen artifact

- Prediction file: `results/v1_3/holdout_predictions_unscored.tsv`
- Rows: {len(out)} across {len(summary)} proteins
- SHA256: `{manifest['hash']}` (computed once and cached in `artifacts/cache/file_manifest.json`)
- Method: official DeepPBS checkout commit `8bfb211dd67f02877841f6f33aa493ddf7daedf9`, bundled five-checkpoint `DeepPBS.txt` ensemble
- Score: `log(P_mutant + 1e-12) - log(P_wt + 1e-12)`; sign unchanged
- Mapping: protein chain A and DNA chain B; position i is PDB/DeepPBS first-strand index i; no register or orientation search
- Coverage: every PDB first-strand position and all three substitutions, independent of assay coverage

## Label-free prediction summary

| protein | sequence | positions | rows | minimum predicted effect | maximum predicted effect |
|---|---|---:|---:|---:|---:|
"""
    for _, row in summary.iterrows():
        report += f"| {row.protein} | `{row.sequence}` | {row.positions} | {row.rows} | {row.prediction_min:.6f} | {row.prediction_max:.6f} |\n"
    report += "\nAll three structures passed chain, duplex, output-shape, probability, and first-strand sequence checks. No experimental label was used for preprocessing, mapping, sign, checkpoint, or failure handling. Evaluation is prohibited until this artifact and report are committed.\n"
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(report, encoding="utf-8")
    print(json.dumps({"status": "frozen_unscored", "rows": len(out), "sha256": manifest["hash"], "summary": summaries}, indent=2))


if __name__ == "__main__":
    main()
