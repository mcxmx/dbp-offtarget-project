from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = ROOT / "artifacts/cache/file_manifest.json"


def record_file(path: Path, source_url: str = "", manifest_path: Path = DEFAULT_MANIFEST) -> dict:
    """Hash a file only when its path, size, or mtime differs from the cache."""
    path = path.resolve()
    stat = path.stat()
    key = path.relative_to(ROOT).as_posix()
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        manifest = {"schema_version": 1, "files": {}}
    cached = manifest["files"].get(key)
    if cached and cached["size"] == stat.st_size and cached["mtime_ns"] == stat.st_mtime_ns:
        return cached

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    entry = {
        "path": key,
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "hash": digest.hexdigest(),
        "hash_algorithm": "sha256",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    if source_url:
        entry["source_url"] = source_url
    manifest["files"][key] = entry
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return entry


def main() -> None:
    parser = argparse.ArgumentParser(description="Update the persistent size/mtime-aware SHA256 manifest.")
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--source-url", default="")
    args = parser.parse_args()
    for path in args.paths:
        entry = record_file(path, args.source_url)
        print(json.dumps(entry, sort_keys=True))


if __name__ == "__main__":
    main()
