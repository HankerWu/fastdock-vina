"""Redock both specified stereoisomers of the reported Ki20227 racemate."""

import argparse
import csv
import itertools
from pathlib import Path
from statistics import median
import sys

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)

from analyze import agg, compare, num, write
from run import ARMS, STATES, SEEDS, run


def main(repo, work, exe):
    sys.path.insert(0, str(repo / "benchmarks/kinase_state_affinity/source"))
    from prepare import ligand

    rows = list(csv.DictReader((work / "panel.csv").open()))
    r = next(r for r in rows if r["target"] == "AURKA" and r["compound"] == "Ki-20227")
    parent = Chem.MolFromSmiles(r["smiles"])
    species = list(
        EnumerateStereoisomers(
            parent, options=StereoEnumerationOptions(onlyUnassigned=True)
        )
    )
    if len(species) != 2:
        raise ValueError("expected one unresolved stereo center")
    tmp = work / "racemate"
    tmp.mkdir(exist_ok=True)
    (tmp / "prepared").mkdir(exist_ok=True)
    panel = []
    for m in species:
        ((_, hand),) = Chem.FindMolChiralCenters(m, includeUnassigned=True)
        row = dict(
            r,
            smiles=Chem.MolToSmiles(m),
            cid="Ki20227" + hand,
            compound="Ki20227" + hand,
        )
        panel.append(row)
        ligand(tmp, row)
    for s in STATES["AURKA"]:
        for ext in (".pdbqt", ".json"):
            (tmp / "prepared" / (s + ext)).write_bytes(
                (work / "prepared" / (s + ext)).read_bytes()
            )
    write(tmp / "panel.csv", panel)
    run(tmp, exe, list(ARMS), 2)
    allcells = list(csv.DictReader((work / "cells.csv").open()))
    mols = sorted({a["cid"] for a in allcells if a["target"] == "AURKA"})
    label = {a["cid"]: float(a["kd_nm"]) for a in rows if a["target"] == "AURKA"}
    results = []
    for arm in ARMS:
        other = {
            (k, s): median(
                float(a["pool_score"])
                for a in allcells
                if a["cid"] == k and a["state"] == s and a["arm"] == arm
            )
            for k in mols
            for s in STATES["AURKA"]
        }
        for s in STATES["AURKA"]:
            e = []
            for row in panel:
                poses = [
                    (tmp / "pools" / f"{s}_{row['cid']}_{seed}_{arm}.pdbqt").read_text()
                    for seed in SEEDS
                ]
                e.append(median(num(p, "REMARK VINA RESULT:") for p in poses))
            mixed = agg(e, "exponential")
            other[r["cid"], s] = mixed
            results.append(
                dict(
                    state=s,
                    arm=arm,
                    species0=panel[0]["cid"],
                    energy0=e[0],
                    species1=panel[1]["cid"],
                    energy1=e[1],
                    mixture=mixed,
                )
            )
        scores = {k: agg([other[k, s] for s in STATES["AURKA"]], "mean") for k in mols}
        pairs = [
            (a, b)
            for a, b in itertools.permutations(mols, 2)
            if label[b] >= 3 * label[a]
        ]
        involved = [(a, b) for a, b in pairs if r["cid"] in (a, b)]
        print(
            arm,
            "all",
            sum(compare(scores[b] - scores[a]) for a, b in pairs),
            len(pairs),
            "involving racemate",
            sum(compare(scores[b] - scores[a]) for a, b in involved),
            len(involved),
        )
    write(tmp / "mixtures.csv", results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for f in ("repo", "work", "exe"):
        p.add_argument(f, type=Path)
    a = p.parse_args()
    main(a.repo.resolve(), a.work.resolve(), a.exe.resolve())
