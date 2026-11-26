"""Fetch the exact external structures used for the water and displacement hypotheses."""

import argparse
import hashlib
import json
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import urlopen


def fetch(out):
    root = Path(__file__).resolve().parent
    out.mkdir(exist_ok=True, parents=True)
    manifest = {}
    for group in (
        "results/donors.json",
        "results/witness_donors.json",
        "holdout_results/donors.json",
    ):
        for accession, sha in json.loads((root / group).read_text()).items():
            if accession in manifest and manifest[accession] != sha:
                raise ValueError(f"Inconsistent digest for {accession}")
            manifest[accession] = sha
    for accession, sha in sorted(manifest.items()):
        target = out / f"{accession}.pdb"
        if target.exists():
            data = target.read_bytes()
        else:
            with urlopen(
                f"https://files.rcsb.org/download/{accession}.pdb", timeout=45
            ) as response:
                data = response.read()
        if hashlib.sha256(data).hexdigest() != sha:
            raise ValueError(f"{accession} differs from the frozen public file")
        if not target.exists():
            with NamedTemporaryFile(dir=out, delete=False) as tmp:
                tmp.write(data)
            Path(tmp.name).replace(target)
    print(f"Verified {len(manifest)} external structures")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    fetch(parser.parse_args().directory)
