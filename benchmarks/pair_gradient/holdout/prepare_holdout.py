"""Prepare prospectively selected androgen receptor cases with the shared pose protocol."""

import argparse
import json
import os
import sys
from pathlib import Path

POSE = Path(__file__).resolve().parents[1] / "pose"
sys.path.insert(0, str(POSE))
import prepare

PRIMARY = {"1E3G": ("AR", "R18", "A"), "2AM9": ("AR", "TES", "A")}
ALTERNATES = {"1XQ3": ("AR", "R18", "A"), "2AMA": ("AR", "DHT", "A")}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    work = parser.parse_args().work
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    prepare.CASES.update(PRIMARY | ALTERNATES)
    ids = list(PRIMARY)
    for pdb in ids:
        prepare.prepare(work, pdb)
    if all(
        json.loads((work / pdb / "preparation.json").read_text())["receptor_status"]
        != "ready"
        for pdb in ids
    ):
        for pdb in ALTERNATES:
            prepare.prepare(work, pdb)
            ids.append(pdb)
    (work / "cases.json").write_text(json.dumps(ids) + "\n")
