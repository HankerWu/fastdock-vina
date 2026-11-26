"""Compare the same top pose before and after two ligand-only geometry relaxations."""

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from posebusters import PoseBusters


def check(case):
    sys.path.insert(0, str(Path(case["repo"]) / "tools"))
    from relax_pose_geometry import first_pose, relax

    mol = first_pose(case["pool"])
    ref = Chem.SDMolSupplier(case["reference"], removeHs=True)[0]
    orig = Chem.RemoveHs(mol)
    pb = PoseBusters(config="dock", max_workers=0)
    rows = []
    for name, k in [("native", None), ("tether10", 10.0), ("free", 0.0)]:
        fail = ""
        try:
            q = mol if k is None else relax(mol, k)
        except ValueError as e:
            fail = str(e)
            q = mol
        heavy = Chem.RemoveHs(q)
        rmsd = rdMolAlign.CalcRMS(heavy, ref, maxMatches=100000)
        df = pb.bust(q, mol_cond=case["protein"])
        checks = dict(df.iloc[0])
        if len(checks) < 15 or any(
            not isinstance(v, (bool, np.bool_)) for v in checks.values()
        ):
            raise ValueError(checks)
        shift = heavy.GetConformer().GetPositions() - orig.GetConformer().GetPositions()
        motion = float(np.sqrt(np.sum(shift**2) / heavy.GetNumAtoms()))
        rows.append(
            dict(
                case=case["id"],
                dataset=case["dataset"],
                query=case["query"],
                arm=name,
                rmsd=rmsd,
                valid=int(all(checks.values())),
                joint=int(all(checks.values()) and rmsd <= 2),
                displacement=motion,
                failed_tests="|".join(k for k, v in checks.items() if not v),
                abstention=fail,
                expected_native=case["expected_native"],
            )
        )
    return rows


def main(repo, pose, hal, output, jobs):
    cases = []

    def cell(dataset, pool, query, ref, protein, rms):
        return dict(
            id=pool.stem,
            dataset=dataset,
            query=query,
            repo=str(repo),
            pool=str(pool),
            reference=str(ref),
            protein=str(protein),
            expected_native=rms,
        )

    for p in sorted(pose.glob("modes-*-250000.csv")):
        for r in csv.DictReader(p.open()):
            if r["rank"] != "1":
                continue
            q, s = r["ligand"], r["receiver"]
            pool = (
                pose
                / "runs_250000"
                / ("_".join((q, s, r["seed"], r["arm"])) + ".pdbqt")
            )
            cases.append(
                cell(
                    "cross",
                    pool,
                    q,
                    pose / "frames" / f"{q}_{s}.sdf",
                    pose / "validity" / f"{s}-condition.pdb",
                    r["rmsd"],
                )
            )
    for r in csv.DictReader((hal / "results/cells.csv").open()):
        d = hal / "data" / r["case"]
        if r["case"].startswith("6YQW_"):
            continue  # BRD9 was used in the cross-docking panel
        cases.append(
            cell(
                "external",
                d / f"poses-{r['seed']}.pdbqt",
                r["case"],
                d / "crystal.sdf",
                d / "polymer.pdb",
                r["native_rmsd"],
            )
        )
    with ProcessPoolExecutor(jobs) as e:
        results = [r for rows in e.map(check, cases) for r in rows]
    for r in results:
        if (
            r["arm"] == "native"
            and abs(r["rmsd"] - float(r["expected_native"])) > 0.015
        ):
            raise ValueError(r)
    with output.open("w") as f:
        w = csv.DictWriter(f, fieldnames=results[0])
        w.writeheader()
        w.writerows(results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("repo", "pose", "hal", "output"):
        p.add_argument(arg, type=Path)
    p.add_argument("--jobs", type=int, default=3)
    a = p.parse_args()
    main(
        a.repo.resolve(), a.pose.resolve(), a.hal.resolve(), a.output.resolve(), a.jobs
    )
