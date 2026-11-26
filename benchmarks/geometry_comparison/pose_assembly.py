"""Paired BCL6 cognate/reciprocal tests of deposited receptor assemblies."""

import argparse
import csv
import itertools
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import gemmi
import numpy as np
from pose_run import HERE, align_with_numbering, load_meta, measure, write_csv

sys.path.insert(0, str(HERE.parents[1] / "tools"))
from build_receptor_assembly import build, coordinates


def mate_transform(deposition, image):
    st = gemmi.read_structure(str(deposition))
    op = gemmi.Op(image["symmetry"])
    orth, frac = (
        np.array(st.cell.orth.mat.tolist()),
        np.array(st.cell.frac.mat.tolist()),
    )
    m = orth @ (np.array(op.rot) / op.DEN) @ frac
    t = orth @ (np.array(op.tran) / op.DEN + image["shift"])
    return m, t


def add_mate(deposition, monomer, receptor, image, output):
    m, t = mate_transform(deposition, image)
    lines = [
        l
        for l in receptor.read_text().splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    ]
    used = {l[21] for l in lines}
    chain = next(i for i in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if i not in used)
    if image["assemblies"]:
        raise ValueError("expected a nonassembly crystal copy")

    for atom in monomer.read_text().splitlines():
        if atom.startswith(("ATOM  ", "HETATM")):
            xyz = m @ coordinates(atom) + t
            lines.append(
                atom[:6]
                + f"{len(lines) + 1:5d}"
                + atom[11:21]
                + chain
                + atom[22:30]
                + "".join(f"{x:8.3f}" for x in xyz)
                + atom[54:]
            )
    output.write_text("\n".join(lines) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--max-evals", type=int, default=120000)
    p.add_argument("--seeds", type=int, nargs="+", default=[63419, 63421, 91571, 91583])
    p.add_argument("--jobs", type=int, default=1)
    args = p.parse_args()
    inputs, out = args.inputs.resolve(), args.output.resolve()
    rows = list(csv.DictReader((HERE / "pose/packing.csv").open()))
    pdbs = sorted(
        r["pdb"]
        for r in rows
        if any(i["assemblies"] for i in json.loads(r["close_images"]))
    )
    meta = {pdb: load_meta(inputs, pdb) for pdb in pdbs}
    policies = {"monomer": [], "dimer": [], "dimer_plus_lattice": []}
    for pdb in pdbs:
        src = inputs / pdb
        deposition = src / (pdb + ".pdb")
        monomer = src / "receptor.pdbqt"
        for pol in policies:
            work = out / pol / pdb
            work.mkdir(parents=True, exist_ok=True)
            (work / "ligand.pdbqt").write_bytes((src / "ligand.pdbqt").read_bytes())
        (out / "monomer" / pdb / "receptor.pdbqt").write_bytes(monomer.read_bytes())
        dimer = out / "dimer" / pdb / "receptor.pdbqt"
        info = build(deposition, monomer, "1", dimer)
        (dimer.parent / "assembly.json").write_text(json.dumps(info, indent=2) + "\n")
        imgs = json.loads(next(r for r in rows if r["pdb"] == pdb)["close_images"])
        nonbio = [i for i in imgs if not i["assemblies"]]
        policies["monomer"].append(pdb)
        policies["dimer"].append(pdb)
        if nonbio:
            if len(nonbio) != 1:
                raise ValueError(
                    "Multiple nonassembly images require explicit treatment"
                )
            add_mate(
                deposition,
                monomer,
                dimer,
                nonbio[0],
                out / "dimer_plus_lattice" / pdb / "receptor.pdbqt",
            )
            policies["dimer_plus_lattice"].append(pdb)
    jobs = []
    for q, r, seed in itertools.product(pdbs, pdbs, args.seeds):
        if meta[q]["target"] != meta[r]["target"]:
            continue
        tr, t, rms, n, shift = (
            align_with_numbering(meta[q], meta[r])
            if q != r
            else (np.eye(3), np.zeros(3), 0.0, 0, 0)
        )
        refs = np.array(meta[q]["mapped"]) @ tr + t
        for pol, receivers in policies.items():
            if r in receivers:
                work = out / pol
                (work / "poses").mkdir(exist_ok=True)
                job = (
                    work,
                    work / "poses",
                    args.binary.resolve(),
                    meta[q],
                    meta[r],
                    refs,
                    rms,
                    n,
                    shift,
                    seed,
                    "armijo",
                    args.max_evals,
                )
                jobs.append((pol, job))

    def one(task):
        policy, job = task
        cell, modes = measure(job)
        return dict(policy=policy, **cell), [dict(policy=policy, **m) for m in modes]

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(one, jobs))
    write_csv(out / "cells.csv", [v for v, _ in results])
    write_csv(out / "modes.csv", [m for _, modes in results for m in modes])


if __name__ == "__main__":
    main()
