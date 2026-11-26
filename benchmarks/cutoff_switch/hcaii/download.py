"""Verify the public coordinate, ligand graph, and assay sources."""

import hashlib
import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def fetch_sources(raw):
    raw.mkdir(parents=True, exist_ok=True)
    sources = json.loads(Path(__file__).with_name("sources.json").read_text())

    def fetch(item):
        name, source = item
        path = raw / name
        if (
            path.exists()
            and hashlib.sha256(path.read_bytes()).hexdigest() == source["sha256"]
        ):
            return
        with urllib.request.urlopen(source["url"], timeout=60) as response:
            blob = response.read()
        if hashlib.sha256(blob).hexdigest() != source["sha256"]:
            raise ValueError(f"Source revision differs from this experiment: {name}")
        path.write_bytes(blob)

    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(fetch, sources.items()))
