"""Find crystal-mate protein atoms touching each experimental ligand."""

import argparse
import csv
import itertools
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree


def coords(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def audit(work, pdb):
    meta = json.loads((work / pdb / "metadata.json").read_text())
    path = work / pdb / (pdb + ".pdb")
    lines = path.read_text().splitlines()
    st = gemmi.read_structure(str(path))
    heavy = [
        l
        for l in lines
        if l.startswith(("ATOM  ", "HETATM"))
        and l[76:78].strip() not in ("H", "D")
        and float(l[54:60]) > 0
    ]
    lig = [
        l
        for l in heavy
        if l.startswith("HETATM")
        and l[17:20].strip() == meta["ligand"]
        and l[21:27] == meta["ligand_instance"]
    ]
    labels = sorted({l[16] for l in lig} - {" "})
    lig = [l for l in lig if l[16] in (" ", labels[0] if labels else " ")]
    protein = [l for l in heavy if l.startswith("ATOM  ")]
    tree = cKDTree([coords(l) for l in lig])
    frac = np.array(st.cell.frac.mat.tolist())
    orth = np.array(st.cell.orth.mat.tolist())
    lfrac = tree.data @ frac.T
    pfrac = np.array([coords(l) for l in protein]) @ frac.T
    radius = 40.0
    span = radius * np.linalg.norm(frac, axis=1)
    closest = (float("inf"),)
    contacts = set()
    images = {}
    for op in st.find_spacegroup().operations():
        symmetry = pfrac @ (np.array(op.rot) / op.DEN).T + np.array(op.tran) / op.DEN
        lo = np.ceil(lfrac.min(0) - span - symmetry.max(0)).astype(int)
        hi = np.floor(lfrac.max(0) + span - symmetry.min(0)).astype(int)
        for shift in itertools.product(*(range(a, b + 1) for a, b in zip(lo, hi))):
            name = op.triplet()
            if name == "x,y,z" and shift == (0, 0, 0):
                continue
            pos = (symmetry + shift) @ orth.T
            distance, ligand_idx = tree.query(pos, distance_upper_bound=radius)
            ix = int(distance.argmin())
            if distance[ix] < closest[0]:
                closest = (
                    float(distance[ix]),
                    name,
                    shift,
                    protein[ix],
                    lig[int(ligand_idx[ix])],
                )
            for j in np.flatnonzero(distance <= 4.0):
                line = protein[j]
                # Multiple occupied altlocs at one site only count once.
                site = (name, shift, line[21:27], line[17:20], line[12:16])
                contacts.add(site)
                key = (name, shift, line[21])
                if key not in images:
                    matrix = orth @ (np.array(op.rot) / op.DEN) @ frac
                    translation = orth @ (np.array(op.tran) / op.DEN + shift)
                    assemblies = []
                    for assembly in st.assemblies:
                        # Compare relative to the selected source-chain identity copy.
                        origins = [
                            gop
                            for g in assembly.generators
                            if meta["chain"] in g.chains
                            for gop in g.operators
                        ]
                        mates = [
                            gop
                            for g in assembly.generators
                            if line[21] in g.chains
                            for gop in g.operators
                        ]
                        for origin, mate in itertools.product(origins, mates):
                            origin_m = np.array(origin.transform.mat.tolist())
                            m = origin_m.T @ np.array(mate.transform.mat.tolist())
                            t = origin_m.T @ (
                                np.array(list(mate.transform.vec))
                                - np.array(list(origin.transform.vec))
                            )
                            if np.allclose(
                                m, matrix, atol=2e-4, rtol=0
                            ) and np.allclose(t, translation, atol=0.02, rtol=0):
                                assemblies.append(assembly.name)
                    images[key] = {
                        "symmetry": name,
                        "shift": list(shift),
                        "chain": line[21],
                        "assemblies": sorted(set(assemblies)),
                        "distance": float("inf"),
                        "sites": set(),
                    }
                img = images[key]
                img["distance"] = min(img["distance"], float(distance[j]))
                img["sites"].add(site)
    if len(closest) == 1:
        raise ValueError(f"No symmetry mate found within {radius} A for {pdb}")
    distance, name, shift, protein_atom, ligand_atom = closest
    nearby = [
        dict(v, sites=len(v["sites"]), shift=[int(i) for i in v["shift"]])
        for v in sorted(images.values(), key=lambda v: (v["distance"], v["chain"]))
    ]
    return {
        "pdb": pdb,
        "target": meta["target"],
        "nearest_mate_A": distance,
        "close_images": json.dumps(nearby, separators=(",", ":")),
        "mate_sites_le4A": len(contacts),
        "mate_symmetry": name,
        "mate_shift": ",".join(map(str, shift)),
        "mate_atom": protein_atom[12:27],
        "ligand_atom": ligand_atom[12:27],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = [
        audit(args.work, c["pdb"])
        for c in json.loads((args.work / "manifest.json").read_text())
        if c["status"] == "ready"
    ]
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
