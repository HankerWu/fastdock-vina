"""Read the crystallographic coordinating ligand atom from each PDB LINK record."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from prep import records, xyz


def main(work):
    meta = json.loads((work / "prepared" / "metadata.json").read_text())
    for directory in [work / "runs", work / "cognate"]:
        out = []
        for row in csv.DictReader((directory / "rows.csv").open()):
            rec, lig, seed, arm = (row[k] for k in ["rec", "lig", "seed", "arm"])
            zinc = xyz(
                next(
                    l
                    for l in records(work / "raw", rec)
                    if l.startswith("HETATM")
                    and l[17:20].strip() == "ZN"
                    and l[21] == "A"
                )
            )
            link = next(
                l
                for l in records(work / "raw", lig)
                if l.startswith("LINK")
                and l[17:20].strip() == "ZN"
                and l[47:50] == meta[lig]["ccd"]
                and int(l[52:56]) == 303
            )
            zn_atom = link[42:46].strip()
            name_index = meta[lig]["names"].index(zn_atom)
            ser = meta[lig]["qt_serial"][name_index]
            reference = np.array(meta[lig]["references"][rec]["coords"])
            box = [
                float(l.split("=")[1])
                for l in (work / "prepared" / f"{rec}_box.txt").read_text().splitlines()
            ]
            halfsize = np.array(box[3:]) / 2
            margin = (halfsize - abs(reference - box[:3])).min()
            poses = (
                (directory / f"{rec}_{lig}_{seed}_{arm}.pdbqt")
                .read_text()
                .split("ENDMDL")
            )
            for num, part in enumerate(poses, 1):
                if not part.strip():
                    continue
                l = next(
                    l
                    for l in part.splitlines()
                    if l.startswith("ATOM") and int(l[6:11]) == ser
                )
                out.append(
                    {
                        "rec": rec,
                        "lig": lig,
                        "seed": seed,
                        "arm": arm,
                        "model": num,
                        "Zn_N_dist": round(np.linalg.norm(xyz(l) - zinc), 5),
                        "aligned_ref_Zn_N_dist": round(
                            np.linalg.norm(reference[name_index] - zinc), 5
                        ),
                        "aligned_ref_box_margin": round(margin, 4),
                    }
                )
        with (directory / "zinc.csv").open("w") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(out)
        for arm in ["native", "smooth"]:
            top = [r for r in out if r["arm"] == arm and r["model"] == 1]
            with_contact = {
                (r["rec"], r["lig"], r["seed"])
                for r in out
                if r["arm"] == arm and r["Zn_N_dist"] <= 2.5
            }
            print(
                directory.name,
                arm,
                "contacts (top,oracle,n)",
                sum(r["Zn_N_dist"] <= 2.5 for r in top),
                len(with_contact),
                len(top),
                "median-top",
                round(np.median([r["Zn_N_dist"] for r in top]), 3),
                "ref-dist range",
                min(r["aligned_ref_Zn_N_dist"] for r in top),
                max(r["aligned_ref_Zn_N_dist"] for r in top),
                "min-margin",
                min(r["aligned_ref_box_margin"] for r in top),
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    main(parser.parse_args().work)
