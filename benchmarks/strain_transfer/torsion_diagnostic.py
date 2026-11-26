"""Post-hoc test of adding only the MMFF bonded torsion difference to native scoring."""

import argparse
import itertools
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from evaluate import SHARED, fixed_energy, panels, physical, write
from rdkit import Chem
from rdkit.Chem import AllChem


def torsion_energy(mol, conf=-1):
    props = physical.mmff_properties(mol)
    for term in ("Angle", "Bond", "Ele", "Oop", "StretchBend", "VdW"):
        getattr(props, "SetMMFF" + term + "Term")(False)
    return AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf).CalcEnergy()


def reference_torsion(mol):
    template = Chem.Mol(mol)
    template.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = physical.GAS_EMBED_SEED
    params.numThreads = 1
    energies = {}
    for conf in AllChem.EmbedMultipleConfs(
        template, numConfs=physical.GAS_CONFORMERS, params=params
    ):
        try:
            energies[conf] = physical.minimized_energy(
                template, physical.mmff_properties(template), conf, ()
            )
        except ValueError:
            continue
    if not energies:
        raise ValueError("No reference conformer")
    return torsion_energy(template, min(energies, key=energies.get))


def process(job):
    (target, lig), cells, lookup = job
    ref = reference_torsion(
        physical.decode_ligand(physical.read_poses(cells[0]["path"])[0])[0]
    )
    rows = []
    for c in cells:
        group = []
        for pose, rmsd in zip(
            physical.read_poses(c["path"]),
            map(float, c["pool_rmsds"].split(",")),
            strict=True,
        ):
            mol, _ = physical.decode_ligand(pose)
            gas = fixed_energy(mol, physical.mmff_properties(mol))
            old = lookup[(target, lig, c["receptor"], c["seed"], str(pose.model))]
            if abs(gas - float(old["eps1_fixed"])) > 1e-7:
                raise ValueError("Underlying MMFF energy changed")
            delta = torsion_energy(mol) - ref
            score = (
                pose.native_score + float(old["scale"]) * physical.STRAIN_WEIGHT * delta
            )
            group.append(
                {
                    "target": target,
                    "receptor": c["receptor"],
                    "seed": c["seed"],
                    "ligand": lig,
                    "model": pose.model,
                    "rmsd": rmsd,
                    "score": score,
                    "torsion_delta": delta,
                    "measurement": float(old["measurement"]),
                }
            )
        rows.append(min(group, key=lambda r: (r["score"], r["model"])))
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path)
    p.add_argument("results", type=Path)
    p.add_argument("--workers", type=int, default=2)
    args = p.parse_args()
    from evaluate import load

    lookup = {
        (r["target"], r["ligand"], r["receptor"], r["seed"], r["model"]): r
        for r in load(args.results / "poses.csv")
    }
    groups = defaultdict(list)
    for key, cell in panels(args.evidence):
        groups[key].append(cell)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        rows = [
            r
            for block in executor.map(
                process, [(key, cells, lookup) for key, cells in groups.items()]
            )
            for r in block
        ]
    pairs = []
    for target, recs in SHARED.items():
        for rec, seed in itertools.product(sorted(recs), ("271828", "314159")):
            group = [
                r
                for r in rows
                if (r["target"], r["seed"], r["receptor"]) == (target, seed, rec)
            ]
            for a, b in itertools.combinations(
                sorted(group, key=lambda r: r["measurement"]), 2
            ):
                if b["measurement"] >= 3 * a["measurement"]:
                    gap = b["score"] - a["score"]
                    pairs.append(
                        {
                            "target": target,
                            "receptor": rec,
                            "seed": seed,
                            "strong": a["ligand"],
                            "weak": b["ligand"],
                            "gap": gap,
                            "credit": float(gap > 0) + 0.5 * (gap == 0),
                        }
                    )
    write(args.results / "torsion_cells.csv", rows)
    write(args.results / "torsion_pairs.csv", pairs)
