from __future__ import annotations

import hashlib
from datetime import date
from pathlib import Path

import esm
import pandas as pd
import torch

from src.utils import ensure_dir, project_root


ROOT = project_root()
SEQUENCE_PATH = ROOT / "metadata" / "v0_3" / "designed_dbp_sequences.csv"
TARGET_PATH = ROOT / "metadata" / "v0_5" / "designed_target_manifest_v0_5.csv"
OUTPUT_DIR = ROOT / "data" / "interim" / "v0_5_local"
METADATA_DIR = ROOT / "metadata" / "v0_5_local"
MODEL_NAME = "esm2_t12_35M_UR50D"
LAYER = 12


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    ensure_dir(OUTPUT_DIR)
    ensure_dir(METADATA_DIR)
    sequences = pd.read_csv(SEQUENCE_PATH).set_index("protein_id")
    targets = pd.read_csv(TARGET_PATH)
    designed_ids = sorted(targets["dbp_id"].astype(str).unique())
    missing = sorted(set(designed_ids) - set(sequences.index.astype(str)))
    if missing:
        raise ValueError(f"Missing designed protein sequences: {missing}")

    audit_rows = []
    batch = []
    for protein_id in designed_ids:
        sequence = str(sequences.loc[protein_id, "protein_sequence"]).strip().upper()
        if not sequence or any(residue not in set("ACDEFGHIKLMNPQRSTVWY") for residue in sequence):
            raise ValueError(f"Invalid designed protein sequence for {protein_id}")
        audit_rows.append(
            {
                "dbp_id": protein_id,
                "sequence_length": len(sequence),
                "sequence_source": str(sequences.loc[protein_id, "source_url"]),
                "sequence_type": "designed_construct",
                "experimental_construct_match": (
                    "reported designed construct sequence; exact PBM expression tag/boundary "
                    "metadata is not separately resolved"
                ),
                "complete_for_model": True,
                "sequence_confidence": str(sequences.loc[protein_id, "sequence_confidence"]),
                "notes": str(sequences.loc[protein_id, "notes"]),
            }
        )
        batch.append((protein_id, sequence))

    pd.DataFrame(audit_rows).to_csv(METADATA_DIR / "designed_protein_sequence_audit.csv", index=False)

    torch.set_num_threads(2)
    model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
    model.eval()
    converter = alphabet.get_batch_converter()
    _, _, tokens = converter(batch)
    with torch.no_grad():
        outputs = model(tokens, repr_layers=[LAYER], return_contacts=False)
    representations = outputs["representations"][LAYER].cpu()

    rows = []
    manifest_rows = []
    for batch_index, (protein_id, sequence) in enumerate(batch):
        # ESM output has BOS at token 0 and EOS immediately after the residue tokens.
        residue_embedding = representations[batch_index, 1 : len(sequence) + 1]
        if residue_embedding.shape != (len(sequence), 480):
            raise ValueError(f"Unexpected residue embedding shape for {protein_id}: {tuple(residue_embedding.shape)}")
        sequence_hash = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        for residue_index, residue_vector in enumerate(residue_embedding.numpy(), start=1):
            row = {
                "protein_id": protein_id,
                "residue_index": residue_index,
                "residue": sequence[residue_index - 1],
            }
            row.update({f"emb_{dimension:03d}": float(value) for dimension, value in enumerate(residue_vector)})
            rows.append(row)
        manifest_rows.append(
            {
                "dbp_id": protein_id,
                "checkpoint_name": MODEL_NAME,
                "repr_layer": LAYER,
                "embedding_dim": 480,
                "raw_sequence_length": len(sequence),
                "embedded_sequence_length": int(residue_embedding.shape[0]),
                "sequence_sha256": sequence_hash,
                "token_handling": "BOS/EOS removed; residue rows retained in input order",
                "protein_lm_frozen": True,
                "local_relative_path": "data/interim/v0_5_local/designed_residue_embeddings_esm2_t12_35M_UR50D.parquet",
                "source_checkpoint": "fair-esm esm2_t12_35M_UR50D",
                "notes": "Frozen residue-level representation for Phase 5A; no fine-tuning.",
            }
        )

    output_path = OUTPUT_DIR / "designed_residue_embeddings_esm2_t12_35M_UR50D.parquet"
    pd.DataFrame(rows).to_parquet(output_path, index=False)
    embedding_hash = sha256_file(output_path)
    for row in manifest_rows:
        row["embedding_file_sha256"] = embedding_hash
        row["created_date"] = date.today().isoformat()
    pd.DataFrame(manifest_rows).to_csv(METADATA_DIR / "residue_embedding_manifest.csv", index=False)
    print(f"wrote {output_path} rows={len(rows)} sha256={embedding_hash}")


if __name__ == "__main__":
    main()
