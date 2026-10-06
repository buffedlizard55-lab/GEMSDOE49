#!/usr/bin/env python3
"""Restore integrity-pinned input rasters from the public owner-hosted mirror.

This does not authenticate the mirrored bytes as organizer originals. The official
DrivenData data page is login-walled; see data/README.md and data/manifest.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gh_raw(repo: str, ref: str, source_path: str, destination: Path) -> None:
    if shutil.which("gh") is None:
        raise RuntimeError("GitHub CLI 'gh' is required to restore the public mirror")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp = destination.with_name(destination.name + ".partial")
    command = [
        "gh", "api", f"repos/{repo}/contents/{source_path}?ref={ref}",
        "-H", "Accept: application/vnd.github.raw",
    ]
    try:
        with temp.open("wb") as output:
            subprocess.run(command, stdout=output, check=True)
        temp.replace(destination)
    except Exception:
        temp.unlink(missing_ok=True)
        raise


def restore(entry: dict, root: Path) -> dict:
    destination = root / entry["dest"]
    if destination.is_file():
        if destination.stat().st_size == entry["bytes"] and sha256(destination) == entry["sha256"]:
            return {"id": entry["id"], "status": "already-verified", "bytes": entry["bytes"], "sha256": entry["sha256"]}
        destination.unlink()

    started = time.monotonic()
    assembling = destination.with_name(destination.name + ".assembling")
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        if "parts" in entry:
            with assembling.open("wb") as combined:
                for index, source_path in enumerate(entry["parts"]):
                    part_path = root / "raw_parts" / f"{entry['id']}-{index:03d}.part"
                    gh_raw(entry["repo"], entry["ref"], source_path, part_path)
                    with part_path.open("rb") as part:
                        shutil.copyfileobj(part, combined, 1 << 20)
                    part_path.unlink(missing_ok=True)
        else:
            gh_raw(entry["repo"], entry["ref"], entry["path"], assembling)

        actual_size = assembling.stat().st_size
        actual_hash = sha256(assembling)
        if actual_size != entry["bytes"] or actual_hash != entry["sha256"]:
            raise RuntimeError(
                f"Integrity check failed for {entry['id']}: expected {entry['bytes']} bytes / "
                f"{entry['sha256']}, got {actual_size} / {actual_hash}"
            )
        assembling.replace(destination)
    except Exception:
        assembling.unlink(missing_ok=True)
        raise
    return {
        "id": entry["id"], "status": "restored-and-hash-verified",
        "bytes": destination.stat().st_size, "sha256": sha256(destination),
        "seconds": round(time.monotonic() - started, 2),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-dir", type=Path, default=ROOT / ".cache" / "gems_data")
    parser.add_argument("--only", default="", help="comma-separated manifest IDs")
    args = parser.parse_args()
    root = args.target_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_text())
    selected = {item for item in args.only.split(",") if item}
    receipt = {
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "target_dir": str(root),
        "mirror_provenance": "public owner-hosted GitHub mirror; hashes authenticate mirror consistency only",
        "files": [],
    }
    for entry in manifest["files"]:
        if selected and entry["id"] not in selected:
            continue
        result = restore(entry, root)
        receipt["files"].append(result)
        print(f"[{result['status']}] {result['id']}: {result['bytes']:,} bytes ({result.get('seconds', 0)} s)", flush=True)
    unknown = selected - {entry["id"] for entry in manifest["files"]}
    if unknown:
        raise SystemExit(f"Unknown manifest IDs: {', '.join(sorted(unknown))}")
    receipt["all_verified"] = all(item["status"] in {"restored-and-hash-verified", "already-verified"} for item in receipt["files"])
    receipt_path = root / "restore_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"\nReceipt: {receipt_path}\nALL_VERIFIED={receipt['all_verified']}")
    return 0 if receipt["all_verified"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as error:
        raise SystemExit(f"GitHub API fetch failed (exit {error.returncode}); no unverified data accepted")
