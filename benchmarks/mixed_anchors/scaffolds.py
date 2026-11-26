"""Run the same native/guided split on the published scaffold-anchor matrices."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from shutil import copyfile
import sys
from dock import ROOT, run, write

sys.path.insert(0, str(ROOT / "benchmarks/mapping_ensemble"))
from score_matrix import evaluate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for n in ("base", "binary", "out"):
        p.add_argument(n, type=Path)
    p.add_argument("--parity", choices=("odd", "even"), default="even")
    p.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    routing = []
    for panel, name in (
        ("brd4", "fdv-brd4"),
        ("t4l", "fdv-pro-fixtures"),
        ("hcaii", "fdv-second"),
    ):
        source = a.base / ("fdv-map-" + panel)
        prepared = a.base / name / "prepared"
        destination = a.out / panel
        destination.mkdir(parents=True, exist_ok=True)
        cases = []
        seeds = set()
        for prior in sorted(source.glob("*_multi.log")):
            rec, lig, seed, _ = prior.stem.split("_")
            seed = int(seed)
            seeds.add(seed)
            native = source / f"{rec}_{lig}_{seed}_native.log"
            site = source / f"{rec}_{lig}.multi"
            for kind in ("single", "multi"):
                copyfile(
                    source / f"{rec}_{lig}.{kind}", destination / f"{rec}_{lig}.{kind}"
                )
            cases.append(
                (
                    panel,
                    rec,
                    lig,
                    seed,
                    prepared / f"{lig}_ligand.pdbqt",
                    prepared / f"{rec}_receptor.pdbqt",
                    prepared / f"{rec}_box.txt",
                    site,
                    native,
                    prior,
                )
            )
        with ThreadPoolExecutor(a.jobs) as workers:
            routing.extend(
                workers.map(
                    lambda c: run(c, a.binary.resolve(), a.out, a.parity), cases
                )
            )
        evaluate(panel, prepared, destination, sorted(seeds), ["mixed"])
        print(panel, len(cases), flush=True)
    write(a.out / "routing_scaffolds.csv", routing)


if __name__ == "__main__":
    main()
