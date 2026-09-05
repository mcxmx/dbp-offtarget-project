from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from src.v0_5_local_training import LocalConfig


ROOT = Path(__file__).resolve().parents[1]


def test_local_config_registers_previously_unseen_fold_and_frozen_protocol():
    config = LocalConfig.from_json()
    assert config.smoke_split_name == "protein_cluster_loco"
    assert config.smoke_fold_id == "protein_cluster_loco_fold_2"
    assert config.seed == 42
    assert config.pair_count_per_protein == 512
    assert config.epochs == 18


def test_residue_embedding_manifest_has_provenance_and_special_token_rule():
    manifest = pd.read_csv(ROOT / "metadata" / "v0_5_local" / "residue_embedding_manifest.csv")
    assert set(manifest["dbp_id"]) == {"DBP1", "DBP3", "DBP5", "DBP6", "DBP9", "DBP35", "DBP48"}
    assert manifest["protein_lm_frozen"].eq(True).all()
    assert manifest["embedding_dim"].eq(480).all()
    assert manifest["token_handling"].str.contains("BOS/EOS removed").all()
    assert manifest["sequence_sha256"].str.len().eq(64).all()
    assert manifest["local_relative_path"].str.contains(r"^data/interim/").all()


def test_residue_cache_has_one_embedding_row_per_input_residue():
    audit = pd.read_csv(ROOT / "metadata" / "v0_5_local" / "designed_protein_sequence_audit.csv")
    cache_path = ROOT / "data" / "interim" / "v0_5_local" / "designed_residue_embeddings_esm2_t12_35M_UR50D.parquet"
    cache = pd.read_parquet(cache_path)
    assert set(cache["protein_id"]) == set(audit["dbp_id"])
    assert cache["residue_index"].ge(1).all()
    assert cache.filter(regex=r"^emb_\d{3}$").shape[1] == 480
    assert cache.groupby("protein_id").size().to_dict() == audit.set_index("dbp_id")["sequence_length"].to_dict()


def test_local_config_does_not_change_frozen_global_config():
    old = json.loads((ROOT / "metadata" / "v0_5" / "v0_5_model_config.json").read_text(encoding="utf-8"))
    local = json.loads((ROOT / "metadata" / "v0_5_local" / "local_model_config.json").read_text(encoding="utf-8"))
    assert old["smoke_fold_id"] == "protein_cluster_loco_fold_1"
    assert local["smoke_fold_id"] == "protein_cluster_loco_fold_2"


def test_frozen_primary_result_hashes_unchanged():
    manifest = (ROOT / "results" / "v0_5" / "PRIMARY_RESULTS_FROZEN_MANIFEST.txt").read_text(encoding="utf-8")
    for line in manifest.splitlines():
        if not line.startswith("results/v0_5/"):
            continue
        relative, expected = line.split(" ", 1)
        digest = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        assert digest == expected
