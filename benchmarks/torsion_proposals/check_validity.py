"""Run PoseBusters docking physical validity checks on top poses in their actual receiver frames."""

import argparse
import csv
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


def test_group(args):
    work, source, receiver, paths = args
    import numpy as np
    from meeko import PDBQTMolecule, RDKitMolCreate
    from posebusters import PoseBusters
    from rdkit import Chem

    pb = PoseBusters(config="dock", max_workers=0)
    # The full deposition contains the ligand that was removed before docking.
    # Use the conditioned protein without crystal ligand, retaining receptor ions.
    conditioned = work / "validity" / (receiver + "-condition.pdb")
    rows = []
    for path in paths:
        mols = RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(str(path), skip_typing=True)
        )
        if len(mols) != 1 or mols[0] is None:
            raise ValueError(path)
        q = Chem.Mol(mols[0])
        first = Chem.Conformer(q.GetConformer())
        q.RemoveAllConformers()
        q.AddConformer(first, assignId=True)
        # All checks in 'dock' use only the predicted ligand and the conditioning protein.
        df = pb.bust(q, mol_cond=str(conditioned), full_report=False)
        tests = {str(k): bool(v) for k, v in df.iloc[0].items()}
        if (
            len(df) != 1
            or len(tests) < 15
            or any(not isinstance(v, (bool, np.bool_)) for v in df.iloc[0].values)
        ):
            raise ValueError(tests)
        rows.append(
            dict(
                cap=path.parent.name.split("_")[1],
                ligand=source,
                receiver=receiver,
                seed=path.stem.split("_")[2],
                arm=path.stem.split("_")[3],
                valid=int(all(tests.values())),
                failed="|".join(k for k, v in tests.items() if not v),
            )
        )
    return rows


def main(work, jobs, bridge, deep):
    meta = json.loads((work / "metadata.json").read_text())
    saved = work / "validity"
    saved.mkdir(exist_ok=True)
    for receptor, m in meta.items():
        lines = (work / receptor / "chain-A.pdb").read_text().splitlines()
        metals = [
            line
            for line in (work / "raw" / (receptor + ".pdb")).read_text().splitlines()
            if line.startswith("HETATM")
            and line[21] == "A"
            and line[17:20].strip() in {"MG", "ZN"}
        ]
        (saved / (receptor + "-condition.pdb")).write_text(
            "\n".join(line for line in lines if not line.startswith("END"))
            + "\n"
            + "\n".join(metals)
            + "\nEND\n"
        )
    groups = {}
    for csvfile in sorted(work.glob("modes-*.csv")):
        if ("-bridge" in csvfile.stem) != bridge or ("-250000" in csvfile.stem) != deep:
            continue
        for r in csv.DictReader(csvfile.open()):
            if r["rank"] != "1":
                continue
            q, t = r["ligand"], r["receiver"]
            name = "_".join([q, t, r["seed"], r["arm"]]) + ".pdbqt"
            groups.setdefault((q, t), []).append(work / ("runs_" + r["cap"]) / name)
    with ProcessPoolExecutor(jobs) as ex:
        rows = [
            r
            for batch in ex.map(
                test_group, [(work, q, t, paths) for (q, t), paths in groups.items()]
            )
            for r in batch
        ]
    with (
        work
        / (
            "physical-validity-deep.csv"
            if deep
            else "physical-validity-bridge.csv"
            if bridge
            else "physical-validity.csv"
        )
    ).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--bridge", action="store_true")
    p.add_argument("--deep", action="store_true")
    a = p.parse_args()
    main(a.work, a.jobs, a.bridge, a.deep)
