"""Prepare, dock and archive the preregistered four receptor/anomer conditions."""

import argparse
import hashlib
import itertools
import json
from importlib.metadata import version
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH / "hydration_transfer"))
from panel import cross, timed_process, write_csv

CASES = {"3IDH": ["GLC"], "1V4S": ["GLC"], "4ZW9": ["GLC", "BGC"]}


def prepare(w):
    raw = w / "raw"
    raw.mkdir(exist_ok=True, parents=True)
    panels = {}
    for pdb, compounds in CASES.items():
        d = w / pdb
        d.mkdir(exist_ok=True)
        path = raw / f"{pdb}.pdb"
        source = f"https://files.rcsb.org/download/{pdb}.pdb"
        checksum = cross.fetch(path, source)
        lines = path.read_text().splitlines()
        for index, comp in enumerate(compounds):
            hetero = [
                l
                for l in lines
                if l.startswith("HETATM") and l[76:78].strip() not in {"H", "D"}
            ]
            deposited = [l for l in hetero if l[21] == "A" and l[17:20].strip() == comp]
            if len(set(l[22:27] for l in deposited)) != 1:
                raise ValueError(f"Expected one deposition site for {comp} in {pdb}")
            refs = cross.dep_records(lines, comp)
            crystal = np.asarray(list(refs[0].values()))
            ccd = raw / f"{comp}.cif"
            ccd_url = f"https://files.rcsb.org/ligands/download/{comp}.cif"
            ccd_sha = cross.fetch(ccd, ccd_url)
            start = time.perf_counter()
            meta = cross.ligand(ccd, refs, d / f"{comp}.pdbqt")
            lig_s = time.perf_counter() - start
            if index == 0:
                start = time.perf_counter()
                preparation = cross.receptor(lines, d, pdb, crystal)
                rec_s = time.perf_counter() - start
            near = {
                name: float(
                    cdist(
                        [cross.coords(l) for l in hetero if l[17:20].strip() == name],
                        crystal,
                    ).min()
                )
                for name in sorted(set(l[17:20].strip() for l in hetero))
            }
            panels[f"{pdb}-{comp}"] = dict(
                pdb=pdb,
                comp=comp,
                center=crystal.mean(axis=0).tolist(),
                **meta,
                **preparation,
                occupancy=sorted(set(float(l[54:60]) for l in deposited)),
                closest_by_hetero_type_A=near,
                sources={source: checksum, ccd_url: ccd_sha},
                ligand_s=lig_s,
                receptor_s=rec_s if index == 0 else 0.0,
            )
    (w / "metadata.json").write_text(json.dumps(panels, indent=2) + "\n")
    (w / "packages.json").write_text(
        json.dumps(
            {
                p: version(p)
                for p in (
                    "numpy",
                    "scipy",
                    "meeko",
                    "rdkit",
                    "openmm",
                    "pdbfixer",
                    "gemmi",
                )
            },
            indent=2,
        )
        + "\n"
    )


def dock(w, bins):
    panels = json.loads((w / "metadata.json").read_text())
    (w / "runs").mkdir(exist_ok=True)
    costs = []
    for name, m in panels.items():
        for limit, seed, arm in itertools.product((3000, 12000), (43103, 43117), bins):
            tag = f"{name}_{limit}_{seed}_{arm}"
            args = [
                "--receptor",
                f"{m['pdb']}/receptor.pdbqt",
                "--ligand",
                f"{m['pdb']}/{m['comp']}.pdbqt",
                "--cpu",
                "1",
                "--seed",
                str(seed),
                "--exhaustiveness",
                "4",
                "--max_evals",
                str(limit),
                "--num_modes",
                "9",
                "--energy_range",
                "100",
                "--out",
                f"runs/{tag}.pdbqt",
            ]
            for axis, value in zip("xyz", m["center"]):
                args.extend(["--center_" + axis, str(value), "--size_" + axis, "24"])
            if arm == "pro":
                args.append("--search_counters")
            log = w / "runs" / (tag + ".log")
            cost = timed_process([str(bins[arm]), *args], w, log)
            trace = log.read_text()
            it = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", trace)
            ex = re.findall(r"explicit_evals=(\d+)", trace)
            costs.append(
                dict(
                    case=name,
                    cap=limit,
                    seed=seed,
                    arm=arm,
                    **cost,
                    args=" ".join(args),
                    search_evals=sum(int(a) + int(b) for a, b in it) if it else "",
                    explicit_evals=sum(map(int, ex)) if ex else "",
                )
            )
    write_csv(w / "compute.csv", costs)
    (w / "binaries.json").write_text(
        json.dumps(
            {
                k: {
                    "version": subprocess.check_output(
                        [str(b), "--version"], text=True
                    ).strip(),
                    "sha256": hashlib.sha256(b.read_bytes()).hexdigest(),
                }
                for k, b in bins.items()
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    sub = p.add_subparsers(dest="action", required=True)
    sub.add_parser("prepare")
    d = sub.add_parser("dock")
    for key in ("upstream", "pro"):
        d.add_argument("--" + key, required=True, type=Path)
    a = p.parse_args()
    if a.action == "prepare":
        prepare(a.work.resolve())
    else:
        dock(
            a.work.resolve(), {k: getattr(a, k).resolve() for k in ("upstream", "pro")}
        )
