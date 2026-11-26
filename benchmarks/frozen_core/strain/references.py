"""Sensitivity to finite unbound references and ligand-only priors in affinity ranking."""

import argparse
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from pose import physical, read, write
from rdkit import Chem
from rdkit.Chem import AllChem

ARMS = ("native", "unit_8", "unit_32", "prepared_prior_8", "prepared_prior_32", "size")


def reference(job):
    compound, path = job
    start, smi = physical.decode_ligand(
        physical.Pose(0, 0, 0, path.read_text(), None, None)
    )
    conf = start.GetConformer()
    fixed = [a.GetIdx() for a in start.GetAtoms() if a.GetAtomicNum() > 1]
    prep = physical.minimized_energy(
        start, physical.mmff_properties(start), conf.GetId(), fixed
    )
    mol = Chem.Mol(start)
    mol.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed, params.numThreads = 49043, 1
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=32, params=params)
    props = physical.mmff_properties(mol)
    values = []
    for i in ids:
        try:
            values.append(physical.minimized_energy(mol, props, i, ()))
        except ValueError:
            continue
    if not values:
        raise ValueError("No alternate unbound conformer converged")
    return {
        "compound": compound,
        "smiles": smi,
        "heavy_atoms": len(fixed),
        "prepared_mmff": prep,
        "free_32": min(values),
        "converged_32": len(values),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("affinity_work", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    base = read(a.affinity_work / "modes.csv")
    jobs = [
        (comp, a.affinity_work / "prepared" / (comp + ".pdbqt"))
        for comp in sorted({r["compound"] for r in base})
    ]
    with ProcessPoolExecutor(a.jobs) as pool:
        refs = {r["compound"]: r for r in pool.map(reference, jobs)}
    scores = defaultdict(lambda: defaultdict(list))
    for r in base:
        comp = r["compound"]
        ref = refs[comp]
        native, inter, scale, E, free = map(
            float,
            (
                r[x]
                for x in [
                    "native",
                    "inter_only",
                    "flexibility_scale",
                    "fixed_mmff",
                    "free_mmff",
                ]
            ),
        )
        ref["free_8"] = free
        alternatives = {
            "native": native,
            "unit_8": float(r["full_unit"]),
            "unit_32": inter + scale * (E - ref["free_32"]),
            "prepared_prior_8": native + scale * (ref["prepared_mmff"] - free),
            "prepared_prior_32": native
            + scale * (ref["prepared_mmff"] - ref["free_32"]),
            "size": -ref["heavy_atoms"],
        }
        for arm, score in alternatives.items():
            scores[r["cap"], r["target"], r["assay"], r["receptor"], comp, arm][
                r["seed"]
            ].append(score)
    pred = {
        k: statistics.median(min(m) for m in by_seed.values())
        for k, by_seed in scores.items()
    }
    if not all(len(byseed) == 2 for byseed in scores.values()):
        raise ValueError("Expected two independent seeds")
    pairs = []
    for r in read(a.affinity_work / "pairs.csv"):
        if r["method"] != "full_unit":
            continue
        group = (r["cap"], r["target"], r["assay"], r["receptor"])
        for arm in ARMS:
            gap = pred[*group, r["weaker"], arm] - pred[*group, r["stronger"], arm]
            credit = 0.5 if abs(gap) < 1e-10 else float(gap > 0)
            if arm == "unit_8" and credit != float(r["credit"]):
                raise ValueError("Reference replay changed primary result")
            pairs.append(
                {
                    "cap": group[0],
                    "target": group[1],
                    "assay": group[2],
                    "receptor": group[3],
                    "stronger": r["stronger"],
                    "weaker": r["weaker"],
                    "arm": arm,
                    "gap": gap,
                    "credit": credit,
                }
            )
    grouped = defaultdict(list)
    for r in pairs:
        grouped[r["cap"], r["target"], r["assay"], r["receptor"], r["arm"]].append(
            r["credit"]
        )
    write(a.affinity_work / "free_references.csv", list(refs.values()))
    write(a.affinity_work / "control_pairs.csv", pairs)
    write(
        a.affinity_work / "control_summary.csv",
        [
            {
                "cap": k[0],
                "target": k[1],
                "assay": k[2],
                "receptor": k[3],
                "arm": k[4],
                "n": len(v),
                "credit": sum(v),
            }
            for k, v in sorted(grouped.items())
        ],
    )


if __name__ == "__main__":
    main()
