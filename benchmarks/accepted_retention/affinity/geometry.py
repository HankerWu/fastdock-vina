"""Compare every hCA II mode against receptor-aligned crystal geometry and Zn ligation."""

import argparse
import json
import re
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from run import read, write


def xyz(line):
    return [float(line[30 + k : 38 + k]) for k in (0, 8, 16)]


def geometry(work, output):
    prep = work / "ca/prepared"
    meta = json.loads((prep / "metadata.json").read_text())
    rows = []
    for run in read(output / "cells.csv"):
        if run["panel"] != "CAII":
            continue
        rec = run["receptor"]
        lig = run["compound"]
        data = meta[lig]
        ref = Chem.SDMolSupplier(str(prep / (lig + "_graph.sdf")), removeHs=True)[0]
        for i, pt in enumerate(data["references"][rec]["coords"]):
            ref.GetConformer().SetAtomPosition(i, pt)
        pr = Chem.Mol(ref)
        crystal = (work / "ca/raw" / (lig + ".pdb")).read_text().splitlines()
        link = next(
            l
            for l in crystal
            if l.startswith("LINK")
            and l[17:20].strip() == "ZN"
            and l[47:50] == data["ccd"]
            and int(l[52:56]) == 303
        )
        idx = data["names"].index(link[42:46].strip())
        receiver = (work / "ca/raw" / (rec + ".pdb")).read_text().splitlines()
        zinc = np.array(
            xyz(
                next(
                    l
                    for l in receiver
                    if l.startswith("HETATM")
                    and l[17:20].strip() == "ZN"
                    and l[21] == "A"
                )
            )
        )
        dest = work / "raw" / (run["key"] + ".pdbqt")
        modes = dest.read_text().split("ENDMDL")
        for rank, chunk in enumerate(modes[:-1], 1):
            coords = {
                int(l[6:11]): xyz(l)
                for l in chunk.splitlines()
                if l.startswith(("ATOM", "HETATM"))
            }
            for i, ser in enumerate(data["qt_serial"]):
                pr.GetConformer().SetAtomPosition(i, coords[ser])
            rms = rdMolAlign.CalcRMS(
                pr, ref, maxMatches=1000000, symmetrizeConjugatedTerminalGroups=True
            )
            rows.append(
                {
                    "receptor": rec,
                    "compound": lig,
                    "method": run["method"],
                    "rank": rank,
                    "score": re.search(
                        r"^REMARK VINA RESULT:\s+(\S+)", chunk, re.MULTILINE
                    )[1],
                    "rmsd_A": round(rms, 6),
                    "zn_n_A": round(
                        float(np.linalg.norm(zinc - coords[data["qt_serial"][idx]])), 6
                    ),
                }
            )
    write(output / "caii_geometry.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    geometry(args.work, args.output)
