"""Prepare the independent cognate structures with the repository's frozen pipeline."""

import argparse
import json
import os
from pathlib import Path
import sys
from time import perf_counter

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "paired_crossdocking"))
import run as preparation

sys.path.insert(0, str(BENCH / "cubic_grid/human_panel"))
from panel import prepare_entry

sys.path.insert(0, str(BENCH / "cubic_grid/fkbp_ad4"))
import run_official

ENTRIES = [
    ("CDK2", "2BTS", "U32"),
    ("BRD4", "3MXF", "JQ1"),
    ("MDM2", "4HG7", "NUT"),
    ("PPARG", "2PRG", "BRL"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("autogrid", type=Path)
    parser.add_argument(
        "--maps-only",
        action="store_true",
        help="Use archived PDBQTs and verify rebuilt map SHA256 digests",
    )
    args = parser.parse_args()
    work = args.work.resolve()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    for entry in ENTRIES:
        t = perf_counter()
        pdb = entry[1]
        data = (
            json.loads((work / pdb / "preparation.json").read_text())
            if args.maps_only
            else prepare_entry(entry, work, preparation)
        )
        if data["status"] != "ready":
            print(pdb, data["status"], data.get("error"))
            continue
        prep_time = perf_counter() - t
        ligand = (work / pdb / "ligand.pdbqt").read_text().splitlines()
        run_official.TYPES = sorted(
            {l[77:].strip() for l in ligand if l.startswith(("ATOM", "HETATM"))}
        )
        t = perf_counter()
        hashes = run_official.prep_grid(work, work, args.autogrid.resolve(), pdb)
        if args.maps_only:
            recorded = json.loads((work / pdb / "map_provenance.json").read_text())[
                "hashes"
            ]
            if hashes != recorded:
                raise ValueError(f"AutoGrid map hashes differ for {pdb}")
            print(pdb, "all map hashes match")
            continue
        (work / pdb / "map_provenance.json").write_text(
            json.dumps(
                dict(
                    hashes=hashes,
                    types=run_official.TYPES,
                    preparation_seconds=prep_time,
                    autogrid_seconds=perf_counter() - t,
                ),
                indent=2,
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
