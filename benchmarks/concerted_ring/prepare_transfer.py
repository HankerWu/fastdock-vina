"""Prepare the registered human PIM1 transfer using the existing shared recipe."""

import argparse
import json
import os
from pathlib import Path
import sys

CASES = {
    "5EOL": ("PIM1", "5QO", "A"),
    "6QBG": ("CatD", "HVT", "A"),
    "6QBH": ("CatD", "S43", "A"),
    "6QCB": ("CatD", "HWE", "A"),
    "3VHA": ("HSP90a", "VHA", "A"),
    "8H6P": ("CDK2", "WZU", "A"),
}
from urllib.request import urlretrieve

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "benchmarks/macrocycle_proposal/external"))
import prepare_panel


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("ids", nargs="*", choices=list(CASES))
    args = p.parse_args()
    work = args.work.resolve()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    prepare_panel.CASES.update(CASES)
    for pdb in args.ids or CASES:
        prepare_panel.run(REPO, work, pdb)
        record = json.loads((work / pdb / "screening.json").read_text())
        if record["status"] != "ready":
            continue
        listing = (work / "raw" / (pdb + ".pdb")).read_text().splitlines()
        _, ccd, chain = CASES[pdb]
        resid = next(
            l[22:26].strip()
            for l in listing
            if l.startswith("HETATM") and l[17:20] == ccd and l[21] == chain
        )
        sdf = work / "raw" / (pdb + "-bound.sdf")
        if not sdf.exists():
            urlretrieve(
                f"https://models.rcsb.org/v1/{pdb.lower()}/ligand?auth_asym_id={chain}&auth_seq_id={resid}&encoding=sdf&copy_all_categories=false",
                sdf,
            )


if __name__ == "__main__":
    main()
