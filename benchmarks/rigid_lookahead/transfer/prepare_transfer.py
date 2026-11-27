"""Prepare the three frozen PDB/CCD pairs, including any preparation exclusions."""

import argparse
import json
import os
import sys
import time
from importlib.metadata import version
from pathlib import Path

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
import prepare

SYSTEMS = {
    "1RBP": ("RBP4", "RTL", "A"),
    "3FMZ": ("RBP4", "2T1", "A"),
    "1D2S": ("SHBG", "DHT", "A"),
}
DIRECTIONS = [(r, q) for r in ("1RBP", "3FMZ") for q in ("1RBP", "3FMZ")] + [("1D2S", "1D2S")]


def main(work):
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    prepare.CASES.update(SYSTEMS)
    statuses = {}
    for code, (target, ccd, chain) in SYSTEMS.items():
        start = time.monotonic()
        prepare.prepare(work, code)
        meta = json.loads((work / code / "preparation.json").read_text())
        statuses[code] = {
            "target": target, "component": ccd, "chain": chain,
            "status": meta["receptor_status"],
            "wall_s": time.monotonic() - start,
            "pdb": f"https://files.rcsb.org/download/{code}.pdb",
            "pdb_sha256": meta["pdb_sha256"],
            "ccd": f"https://files.rcsb.org/ligands/download/{ccd}.cif",
            "ccd_sha256": meta["ccd_sha256"],
        }
    (work / "sources.json").write_text(json.dumps({
        "systems": statuses,
        "python": sys.version.split()[0],
        "versions": {p: version(p) for p in ("rdkit", "meeko", "pdbfixer", "OpenMM", "gemmi", "numpy", "scipy")},
    }, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    main(parser.parse_args().work.resolve())
