"""Fetch primary assay, coordinate and molecular graph sources and verify SHA256."""

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path


def fetch(work):
    raw = work / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(
        (Path(__file__).resolve().parent / "sources.json").read_text()
    )
    for name, entry in manifest.items():
        path = raw / name
        data = (
            path.read_bytes()
            if path.exists()
            else urllib.request.urlopen(entry["url"], timeout=60).read()
        )
        if hashlib.sha256(data).hexdigest() != entry["sha256"]:
            raise ValueError(f"Source SHA256 mismatch for {name}")
        if not path.exists():
            path.write_bytes(data)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    fetch(p.parse_args().work)
