"""Diagnose within-ligand ordering when a crystal-initialized hCAII pose is explicitly added."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from evaluate import ARMS, TABLES, load, physical, pose_energies, write


def compare(job):
    ligand, binary, evidence = job
    data = evidence / "hcaii"
    prepared = data / "prepared"
    base = prepared / f"{ligand}_ligand.pdbqt"
    receptor = prepared / f"{ligand}_receptor.pdbqt"
    saved = data / "crystal_start"
    bonds = saved / f"{ligand}.bonds"
    historical = [
        r
        for r in load(TABLES / "hcaii/crystal_start.csv")
        if r["ligand"] == ligand and r["arm"] == "native"
    ]

    def scores(path, model):
        # Prepared topology is always reused for Vina; other rankers use the same chemical graph.
        terms = physical.evaluate_coordinates(
            binary, path, receptor, None, 60, prepared_ligand=base, ligand_bonds=bonds
        )
        if model:
            pose = single_saved_pose(path, model)
        else:
            pose = physical.Pose(0, 0, 0, path.read_text(), None, None)
        molecule, _ = physical.decode_ligand(pose)
        gas, screened, basin = pose_energies(molecule)
        inter = terms["inter"]
        w = physical.STRAIN_WEIGHT
        return {
            "native": inter + terms["intra"],
            "inter_only": inter,
            "gas_global": inter + w * gas,
            "screened_global": inter + w * screened,
            "screened_local": inter + w * (screened - basin),
        }

    crystal = scores(saved / f"{ligand}_native_refined.pdbqt", 0)
    comparison = []
    for r in historical:
        seed = int(r["seed"])
        decoys = [
            scores(saved / f"{ligand}_{seed}_{m}.pdbqt", m)
            for m in range(1, int(r["n_modes"]) + 1)
        ]
        for arm in ARMS:
            energies = [p[arm] for p in decoys]
            gap = crystal[arm] - min(energies)
            if arm == "native" and abs(gap - float(r["refined_minus_best"])) > 2e-6:
                raise ValueError(f"Historical crystal gap changed: {r}")
            comparison.append(
                {
                    "ligand": ligand,
                    "seed": seed,
                    "arm": arm,
                    "models": len(decoys),
                    "refined_rmsd": r["refined_RMSD"],
                    "best_decoy_model": energies.index(min(energies)) + 1,
                    "crystal_rank": 1
                    + sum(energy < crystal[arm] for energy in energies),
                    "crystal_minus_best_decoy": gap,
                }
            )
    return comparison


def single_saved_pose(path, model):
    # The historical one-model file preserves its original MODEL id (2..9).
    text = path.read_text().splitlines()
    header, *body, footer = text
    if header.split() != ["MODEL", str(model)] or footer != "ENDMDL":
        raise ValueError(f"Invalid single-pose envelope: {path}")
    body = [
        line
        for line in body
        if not line.startswith(
            ("REMARK VINA ", "REMARK INTER", "REMARK INTRA", "REMARK UNBOUND")
        )
    ]
    return physical.Pose(model, 0, model, "\n".join(body) + "\n", None, None)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--workers", type=int, default=2)
    a = p.parse_args()
    jobs = [
        (lig, a.binary.resolve(), a.evidence)
        for lig in ("6WQ4", "6WQ5", "6WQ7", "6WQ8", "6WQ9")
    ]
    with ProcessPoolExecutor(max_workers=a.workers) as executor:
        rows = [r for group in executor.map(compare, jobs) for r in group]
    a.output.parent.mkdir(parents=True, exist_ok=True)
    write(a.output, rows)


if __name__ == "__main__":
    main()
