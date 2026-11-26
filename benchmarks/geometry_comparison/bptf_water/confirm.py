"""Separate water-site filtering from coordinate relocation on fresh seeds and ligands."""

import argparse
import json
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from registration import match_observed
from water import (
    ANCHOR,
    fixed_water,
    prepare,
    align_with_numbering,
    load_meta,
    measure,
    write_csv,
)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pose_prepare import prepare as prepare_case

EXTERNAL = {"7RWP": "7WZ", "7RWQ": "7XE"}
ARMS = ("dry", "polar", "supported", "registered")
SEEDS = (95701, 95713, 95717)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fixtures", type=Path)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    work = a.work.resolve()
    prepare(a.fixtures, work)
    extra = work / "external"
    extra.mkdir(exist_ok=True)
    queries = {pdb: load_meta(a.fixtures, pdb) for pdb in ANCHOR}
    candidates = [
        dict(pdb=pdb, ligand=ligand, target="BPTF", cohort="external")
        for pdb, ligand in EXTERNAL.items()
    ]
    manifests = [prepare_case((extra, case)) for case in candidates]
    (work / "external_manifest.json").write_text(json.dumps(manifests, indent=2) + "\n")
    for case in manifests:
        if case["status"] == "ready":
            pdb = case["pdb"]
            queries[pdb] = load_meta(extra, pdb)
    registered = {}
    evidence = json.loads((work / "evidence.json").read_text())
    for r in ANCHOR:
        prepared = a.fixtures / r
        sites = match_observed(
            evidence[r]["sites"], prepared / (r + ".pdb"), prepared / "receptor.pdbqt"
        )
        registered[r] = sites
        for arm in ("supported", "registered"):
            out = work / arm / r
            out.mkdir(parents=True, exist_ok=True)
            added = []
            for i, s in enumerate(sites):
                xyz = s["observed_xyz"] if arm == "registered" else s["xyz"]
                added.extend(fixed_water(dict(xyz=xyz), 9001 + 3 * i, 800 + i, True))
            (out / "receptor.pdbqt").write_text(
                (prepared / "receptor.pdbqt").read_text().rstrip()
                + "\n"
                + "".join(l + "\n" for l in added)
            )
    (work / "registered.json").write_text(json.dumps(registered, indent=2) + "\n")
    jobs = []
    for q, qm in queries.items():
        for arm in ARMS:
            directory = work / arm / q
            directory.mkdir(exist_ok=True)
            source = (a.fixtures if q in ANCHOR else extra) / q / "ligand.pdbqt"
            (directory / "ligand.pdbqt").write_bytes(source.read_bytes())
        for r in ANCHOR:
            rm = queries[r]
            rot, tr, rms, n, offset = np.eye(3), np.zeros(3), 0.0, 0, 0
            if q != r:
                rot, tr, rms, n, offset = align_with_numbering(qm, rm)
            ref = np.asarray(qm["mapped"]) @ rot + tr
            for i, seed in enumerate(SEEDS):
                for arm in ARMS[i:] + ARMS[:i]:
                    out = work / arm / "confirm"
                    out.mkdir(exist_ok=True)
                    jobs.append(
                        (
                            work / arm,
                            out,
                            a.binary.resolve(),
                            qm,
                            rm,
                            ref,
                            rms,
                            n,
                            offset,
                            seed,
                            arm,
                            120000,
                        )
                    )
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        results = list(pool.map(measure, jobs))
    write_csv(work / "confirm_cells.csv", [row for row, _ in results])
    write_csv(work / "confirm_modes.csv", [m for _, poses in results for m in poses])


if __name__ == "__main__":
    main()
